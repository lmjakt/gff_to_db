#include "gff_db.h"
#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <set>
#include <map>
#include <stdio.h> // can I use C FILE*

using namespace std;

map<string, int> f_ranks = {
  { "cds", 1 },
  { "exon", 2 },
  { "mrna", 3 },
  { "lnc_rna", 3 },
  { "transcript", 3 },
  { "rna", 3 },
  { "gene", 4 },
  { "region", 5}
};

map<int, string> file_suffix = {
  { 1, "_cds.tsv" },
  { 2, "_exon.tsv" },
  { 3, "_rna.tsv" },
  { 4, "_gene.tsv"},
  { 5, "_region.tsv"}
};

map<int, string> annotation_fields = {
  { 1, "name" },
  { 2, "description" },
  { 3, "product", }
};


void print_feature(FILE* out, const gff_feature& f, size_t region_id){
  // for testing, first do:
  // db_id, parent_id, region_id, feature, id, seqname, parent_name, source, start, end, strand, frame
  //  cout << "print_feature: " << f.feature << " : " << f.id << " : " << f.db_id << " : " << f.parent << " : " << f.parent_id << endl;
  fprintf(out, "%ld %ld %ld %s %s %s %s %s %ld %ld %c %d\n",
	 f.db_id, f.parent_id, region_id, f.feature.c_str(), f.id.c_str(), f.seqname.c_str(), f.parent.c_str(),
	 f.source.c_str(), f.f_start, f.f_end, f.strand, f.frame);
}

void print_annotation(FILE* out, const vector<gff_feature>& gff, int feature_rank){
  for(auto it=gff.begin(); it != gff.end(); ++it){
    if(it->feature_rank == feature_rank && it->annotation.size()){
      for(auto it2=it->annotation.begin(); it2 != it->annotation.end(); ++it2){
	for( auto it3=it2->second.begin(); it3 != it2->second.end(); ++it3){
	  fprintf(out, "%ld\t%d\t%s\n", it->db_id, it2->first, (*it3).c_str());
	}
      }
    }
  }
}

int main(int argc, char **argv){
  if(argc < 3){
    cerr << "Please specify the input gff file and output prefix\n";
    exit(1);
  }
  string o_prefix(argv[2]);
  
  ifstream in(argv[1]);
  if(!in){
    cerr << "Unable to open: " << argv[1] << "\n";
    exit(1);
  }
  // Using a set is about faster as we don't need to sort it afterwards
  set<gff_feature> features;
  string line;
  size_t n_lines = 0;
  while( getline(in, line) ){
    gff_feature feat(line, f_ranks);
    n_lines++;
    if(feat.valid){
      features.insert( feat );
    }
  }
  cerr << "Obtained " << features.size() << " features" << endl
       << "From total of: " << n_lines << endl;
  ///// WARNING /// BUG ALERT ///
  /// This approach will fail if two parents have the same id; 
  /// This means that I will need to group parents by level
  /// as well; For this I should use a specialised struct and
  /// have some form of lookup table;
  // there are three levels that can be parents; rna, gene, region
  gff_parent_collection current_parents( 3 );
  //  map<string, gff_feature> current_parents;
  map<int, size_t> db_ids;
  for(auto it=f_ranks.begin(); it != f_ranks.end(); ++it)
    db_ids[ it->second ] = 0;

  // create the output files:
  // c++ does not let me put ofstreams into a map. Or at least
  // I can't work out how to do it.
  map<int, FILE* > streams;
  string ofile;
  for(auto it=file_suffix.begin(); it != file_suffix.end(); it++){
    ofile = o_prefix + it->second;
    streams[ it->first ] = fopen(ofile.c_str(), "wb");
  }
  // I will also need two files for transcript_exon, and transcript_cds
  string tr_exon_file = o_prefix + "_tr_exon.tsv";
  string tr_cds_file = o_prefix + "_tr_cds.tsv";
  FILE* tr_exon = fopen(tr_exon_file.c_str(), "wb");
  FILE* tr_cds = fopen(tr_cds_file.c_str(), "wb");
  string tr_ann_file = o_prefix + "_tr_annotation.tsv";
  FILE* tr_annotation = fopen(tr_ann_file.c_str(), "wb");

  gff_feature const *last_feature = 0;
  for(set<gff_feature>::iterator it=features.begin(); it != features.end(); ++it){
    // implicit copy needed as we cannot modify *it
    gff_feature feat = gff_feature(*it);
    string feature_id = feat.id; 
    bool repeated = false;
    if(feat.feature_rank < 3 && last_feature) 
      repeated = feat.range_feature_identical(*last_feature);
    last_feature = &(*it);


    vector<gff_feature> discarded = current_parents.prune_parents(feat);
    print_annotation(tr_annotation, discarded, f_ranks["rna"]);
    size_t region_id = 0;
    if(feat.feature == "region"){
      db_ids[ f_ranks["region"] ]++;
      region_id = db_ids[ f_ranks["region"] ];
      feat.db_id = region_id;
      current_parents.insert_parent( feat );
      print_feature( streams[ feat.feature_rank ], feat, region_id );
      if(feat.feature == "region")
	continue;
    }

    // We may still need to define a region id if the the file does not contain regions.
    // find a region identifier; if not found, then make a fake region
    region_id = current_parents.get_parent_db_id( feat.seqname, f_ranks["region"] );
    if(!region_id){
      db_ids[ f_ranks["region"] ]++;
      region_id = db_ids[ f_ranks["region"] ];
      current_parents.insert_parent( gff_feature(region_id, feat.seqname, "region", f_ranks["region"]) );
    }

    // if level is three or lower, a parent should have been defined:
    // if none exists, then we do not want to increment our db_id counter
    // or print anything out.
    gff_feature* parent = current_parents.get_parent( feat );
    if(feat.feature_rank <= 3 && parent == 0){
      cerr << "No parent found for feature: " << feat.id << endl << "\t" << feat.attribute << endl;
      continue;
    }
    feat.parent_id = parent ? parent->db_id : 0;
    // if the feature is a cds or an exon, but the parent is a gene, then create
    // a transcript; print out the information about the transcript as well:
    if(parent != 0){
      if(feat.feature_rank < 3 && parent->feature_rank == f_ranks["gene"]){
	db_ids[ f_ranks["rna"] ]++;
	gff_feature rna(db_ids[f_ranks["rna"]], "rna", f_ranks["rna"], feat);
	print_feature( streams[rna.feature_rank], rna, region_id );
	current_parents.insert_parent(rna); 
	// exons and cds do not actually have single parents, so, this is a bit wrong.
	feat.parent_id = db_ids[ f_ranks["rna"] ];
      }
      int inc = feat.strand == '+' ? 1 : -1;
      parent->child_count += inc;
      if(feat.feature == "exon")
	parent->exon_count += inc;
      if(feat.feature == "cds")
	parent->cds_count += inc;
    }
    // increment the the dabase id:
    // it may be that we should allow rna types to be repeated as they can have 
    // the same range 
    if(!repeated || feat.feature_rank == 3){
      // if the level is three or lower  
      db_ids[ feat.feature_rank ]++;
      feat.db_id = db_ids[ feat.feature_rank ];
      print_feature(streams[feat.feature_rank], feat, region_id);
      if(feat.feature_rank >= 3){
	// gff_parent_collection will convert the id for a region to the region name if this is necessary:
	current_parents.insert_parent(feat);
      }
    }else{ // set the db_id to the current one;
      feat.db_id = db_ids[ feat.feature_rank ];
    }
    if(feat.feature_rank < 3){
      FILE* tr_map = (feat.feature == "cds") ? tr_cds : tr_exon;
      int count = (feat.feature == "cds") ? parent->cds_count : parent->exon_count;
      fprintf(tr_map, "%ld\t%ld\t%d\t%s\n", parent->db_id, feat.db_id, count, feat.id.c_str());
      parent->add_annotation( &feat, annotation_fields );
    }
    if(feat.feature_rank == 3)
      feat.add_annotation( &feat, annotation_fields );
  }
  // and print the last of the annotation:
  if(current_parents.parents[0].size()){
    vector<gff_feature> gff;
    for(auto it=current_parents.parents[0].begin(); it != current_parents.parents[0].end(); ++it)
      gff.push_back( it->second );
    print_annotation( tr_annotation, gff, f_ranks["rna"] );
  }
}

