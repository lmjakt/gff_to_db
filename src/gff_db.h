#ifndef GFF_DB_H
#define GFF_DB_H

#include <string>
#include <map>
#include <set>
#include <vector>
#include <iostream>

#define N_FIELDS 9

// holds an entry from a gff file;
// The fields of a gff are: (taken from (https://www.ensembl.org/info/website/upload/gff.html)
// 
// 1. seqname - name of the chromosome or scaffold; chromosome names can be given with or without the 'chr' prefix. 
// 2. source - name of the program that generated this feature, or the data source (database or project name)
// 3. feature - feature type name, e.g. Gene, Variation, Similarity
// 4. start - Start position* of the feature, with sequence numbering starting at 1.
// 5. end - End position* of the feature, with sequence numbering starting at 1.
// 6. score - A floating point value.
// 7. strand - defined as + (forward) or - (reverse).
// 8. frame - One of '0', '1' or '2'. '0' indicates that the first base of the feature is the first base of a codon, '1' that the second base is the first base of a codon, and so on..
// 9. attribute - A semicolon-separated list of tag-value pairs, providing additional information about each feature.


struct gff_feature {
  bool valid;
  size_t db_id;
  std::string id; // (taken from the attributes)
  std::string parent;
  size_t parent_id;
  // the following are taken directly from the gff:
  std::string seqname;
  std::string source;
  std::string feature;
  // these should never be negative, but don't assume
  unsigned long f_start;
  unsigned long f_end;
  double score;  // usually .
  char strand; // +, -, (.?)
  int frame;
  std::string attribute;
  // split by ; and =
  std::map<std::string, std::string> attributes;

  // This is inferred:
  int feature_rank; // defined as 0, if unknown feature type.
  
  // this is set for rna type (feature_rank == 3) in order to keep track of the
  // child, exons and cds counts
  int child_count;
  int exon_count;
  int cds_count;

  // This field may be set externally to allow the accumulation of annotation;
  std::map< int, std::set<std::string> > annotation;

  // construct from on line of a file:
  // and a map that defines an integer for sorting by feature type
  // So that features with identical coordinates will sort in a deterministic
  // and reasonable order (eg. exon before CDS)
  gff_feature(const std::string& line, const std::map<std::string, int>& feature_ranks);
  gff_feature(size_t db_id, std::string feature, int feature_rank, const gff_feature& src );
  // the following is used to make a feature for a region if these are not defined.
  gff_feature(size_t db_id, std::string seqname, std::string feature, int feature_rank);
  // empty desctructor
  ~gff_feature(){};

  // add annotation:
  void add_annotation(const gff_feature* feat, const std::map<int, std::string>& fields);


  // To allow sorting ensuring that parents are defined before children:
  // This seems like a very inefficient way to compare, but I've not come
  // up with a better one:
  // using std::tie in an external function might be better
  // see: https://stackoverflow.com/questions/3882467/defining-operator-for-a-struct
  //    but we may need to compare by map elements that may not be defined.. 
  inline bool operator< (const gff_feature& other) const
  {
    if(valid == false || other.valid == false)
      return(false);
    if(seqname != other.seqname)
      return(seqname < other.seqname);
    if(f_start != other.f_start)
      return(f_start < other.f_start);
    if(f_end != other.f_end)        // invert the comparison here!
      return(f_end > other.f_end);
    if(feature_rank != other.feature_rank)
      return(feature_rank > other.feature_rank);
    if(strand != other.strand)
      return(strand < other.strand);
    // In the case of a transcript the coordinates may be identical,
    // (since only start and end are given), but the identity may be different
    // In the case of an exon, we may have multiple exons with the same coordinates,
    // same id, but different parents.
    auto t_it = attributes.find("id");
    auto o_it = other.attributes.find("id");
    // Note that valid == true implies that id has been set
    // so we should not need to check this here:
    if(t_it != attributes.end() && o_it != other.attributes.end() && t_it->second != o_it->second)
      return(t_it->second < o_it->second);
    // then compare by parent; this may not be set
    t_it = attributes.find("parent");
    o_it = other.attributes.find("parent");
    if(t_it != attributes.end() && o_it != other.attributes.end() && t_it->second != o_it->second)
      return(t_it->second < o_it->second);
    // return false if nothing else is available
    return(false);
  }

  // Two function that defines whether or not ranges are contained:
  // does "this" contain the other
  // if this has negative start, then, it is considered to have infinite range.
  bool contains(const gff_feature& other) const
  {
    return(other.seqname == seqname && 
	   (f_start < 0 || (other.f_start >= f_start && other.f_end <= f_end)));
  }
  // is "this" contained by the other
  bool contained(const gff_feature& other) const
  {
    return(other.seqname == seqname && f_start >= other.f_start && f_end <= other.f_end);
  }
  // do the two features overlap ?
  bool overlaps(const gff_feature& other) const
  {
    return(other.seqname == seqname && other.f_start <= f_end && other.f_end >= f_start);
  }
  bool range_feature_identical(const gff_feature& other) const
  {
    return(other.seqname == seqname && 
	   other.f_start == f_start && 
	   other.f_end == f_end && 
	   other.feature == feature &&
	   other.strand == strand);
  }
};

inline bool cmp_gff(const gff_feature& a, const gff_feature& b)
{
    if(a.seqname != b.seqname)
      return(a.seqname < b.seqname);
    if(a.f_start != b.f_start)
      return(a.f_start < b.f_start);
    if(a.f_end != b.f_end)        // invert the comparison here!
      return(a.f_end > b.f_end);
    return(a.feature_rank > b.feature_rank);
}

inline std::string str_to_lower(std::string str){
  for(size_t i=0; i < str.size(); ++i){
    if(str[i] >= 'A' && str[i] <= 'Z')
      str[i] = str[i] | 0x20;
  }
  return(str);
}

struct gff_parent_collection {
  // the intended design is that:
  // feature_rank - 3 = vector index
  // but this should be set as part of the constructor;
  std::vector<std::map<std::string, gff_feature> > parents;
  int v_offset;

  gff_parent_collection(int v_offset);
  
  std::vector<gff_feature> prune_parents( const gff_feature& feature );
  void insert_parent( gff_feature feature );
  // Returns a pointer or 0 to the parent; note that this 
  // pointer is temprorary and _must_ not be used after any other
  gff_feature* get_parent(const gff_feature& child);
  size_t get_parent_db_id(std::string parent_id, int rank);
};

#endif
