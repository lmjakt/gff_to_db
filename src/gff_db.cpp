#include "gff_db.h"
#include <fstream>
#include <iostream>
#include <string>
#include <map>
#include <set>
#include <vector>

using namespace std;

map<string, string> split_attributes(const string& line){
  // split by ';', then by '=' to create the attributes
  map<string, string> attr;
  size_t beg = 0;
  size_t end = line.find_first_of(';', beg);
  while(true){
    end = (end == string::npos) ? line.size() : end;
    size_t eq = line.find_first_of('=', beg);
    if(eq != string::npos && eq < end){
      attr[ str_to_lower(line.substr(beg, eq-beg)) ] = 
	str_to_lower(line.substr(eq+1, end-(eq+1)));
    }
    if(end == line.size())
      break;
    beg = end + 1;
    end = line.find_first_of(';', beg);
  }
  return(attr);
}


gff_feature::gff_feature(const string& line, const map<string, int>& feature_ranks){
  // default values:
  db_id = parent_id = 0;
  id = parent = seqname = source = feature = "undefined";
  f_start = f_end = 0;
  exon_count = cds_count = child_count = 0;
  score = 0;
  strand = '.';
  frame = -1;
  // there are meant to be 9 fields (see header)
  vector<string> fields;
  fields.reserve(N_FIELDS);
  // there may be comment lines;
  valid = false;
  if(line.size() < N_FIELDS * 2 || line[0] == '#')
    return;
  size_t beg = 0;
  size_t end = line.find_first_of('\t');
  while(fields.size() < N_FIELDS && end != string::npos){
    fields.push_back( line.substr(beg, end - beg) );
    beg = end + 1;
    end = line.find_first_of('\t', beg);
  }
  if(beg != line.size() && fields.size() < N_FIELDS)
    fields.push_back( line.substr(beg, line.size() - beg));
  if(fields.size() != N_FIELDS)
    return;
  // some of these can throw exceptions; we should consider how
  // to handle these later. For now, it is reasonable to just crash.
  seqname = fields[0];
  source = fields[1];
  feature = str_to_lower(fields[2]);
  f_start = stoul(fields[3]);
  f_end = stoul(fields[4]);
  score = (fields[5] == ".") ? 0 : stod(fields[5]);
  strand = fields[6][0];
  frame = (fields[7] == ".") ? -1 : stod(fields[7]);
  attribute = fields[8];
  attributes = split_attributes(attribute);
  feature_rank = 0;
  // if attributes->id is not set, then valid=false
  if(!attributes.count("id")){
    cerr << "no id attribute set for: " << endl << attribute << endl;
    return;
  }
  //  const map<string, int>::iterator it = feature_ranks.find(feature);
  // use of auto as I don't seem to understand the type that should be
  // used here.
  auto it = feature_ranks.find(feature);
  if(it != feature_ranks.end()){
    valid=true;
    feature_rank = it->second;
  }
  // if the feature is a region, then id will be the region name;
  id = (feature == "region") ? seqname : attributes["id"];
  if(attributes.count("parent"))
    parent = attributes["parent"];
  return;
}

gff_feature::gff_feature(size_t db_id, std::string feature, int feature_rank, const gff_feature& src) :
  valid(true), db_id(db_id), parent(src.parent), seqname(src.seqname), source(src.source), feature(feature),
  f_start(src.f_start), f_end(src.f_end), score(src.score),
  strand(src.strand), frame(src.frame), attributes(src.attributes), feature_rank(feature_rank)
{
  attribute = "";
  exon_count = cds_count =child_count = 0;
}

gff_feature::gff_feature(size_t db_id, std::string seqname, std::string feature, int feature_rank) :
  valid(true), db_id(db_id), id(seqname), seqname(seqname), source("."), feature(feature),
  f_start(0), f_end(-1), score(0), strand('+'), frame(-1), attribute(""), feature_rank(feature_rank)
{
  parent_id=0;
  exon_count = cds_count = child_count = 0;
}

void gff_feature::add_annotation(const gff_feature* feat, const map<int, string>& fields)
{
  for(auto it=fields.begin(); it != fields.end(); ++it){
    auto a_it = feat->attributes.find( it->second );
    if(a_it != feat->attributes.end())
      annotation[ it->first ].insert( a_it->second );
  }
}

gff_parent_collection::gff_parent_collection(int v_offset) :
  v_offset(v_offset)
{
  parents.resize(v_offset);
}

vector<gff_feature> gff_parent_collection::prune_parents( const gff_feature& feature )
{
  vector<gff_feature> discarded;
  for(auto it1=parents.begin(); it1 != parents.end(); ++it1){
    for(auto it2=it1->begin(); it2 != it1->end(); ){
      if(! it2->second.overlaps( feature )){
	discarded.push_back(it2->second);
	it2 = it1->erase(it2);
      }else{
	++it2;
      }
    }
  }
  return(discarded);
}

void gff_parent_collection::insert_parent( gff_feature feature )
{
  int i = feature.feature_rank - v_offset;
  if(i < 0 || i >= (ssize_t)parents.size())
    return;
  // BUG potential; use of literal integer here should be replaced
  string id = (feature.feature == "region") ? feature.seqname : feature.id;
  parents[(unsigned int)i].insert(make_pair(id, feature));
}

gff_feature* gff_parent_collection::get_parent(const gff_feature& child)
{
  // start from a rank one higher than the child to avoid circular
  // links where the id is the same;
  int beg = 1 + child.feature_rank - v_offset;
  beg = beg < 0 ? 0 : beg;
  for(size_t i=(size_t)beg; i < parents.size(); ++i){
    auto it2 = parents[i].find(child.parent);
    if(it2 != parents[i].end())
      return( &(it2->second) );
  }
  return(0);
}

size_t gff_parent_collection::get_parent_db_id(std::string parent_id, int rank){
  int i=rank - v_offset;
  if(i < 0 || i >= (ssize_t)parents.size())
    return(0);
  auto it = parents[i].find(parent_id);
  if(it == parents[i].end())
    return(0);
  return( it->second.db_id );
}
