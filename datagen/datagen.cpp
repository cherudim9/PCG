#include <iostream>
#include <cstdio>
#include <algorithm>
#include <vector>
#include <string>
#include <map>
#include <set>
#include <fstream>
using namespace std;

set<int> net_set;
set<pair<pair<int,int>, int>> node_set;
map<pair<pair<int,int>, int>, int> node_id;
map<int,int> final_id;
map<int,int> father;//union-set
int node_cnt=0;

const double eps=1e-6;

void string_regularize(string &s, bool number_only = false)
{
  while(s.length()>0 && 
	(s[0]==' ' || s[0]=='\n' || s[0]=='\r'
	 || (number_only && (s[0]<'0' || s[0]>'9'))
	 )
	)
    s.erase(0,1);
}

int string_get_number(string &s)
{
  int ret=0;
  while(s.length()>0 && s[0]>='0' && s[0]<='9'){
    ret=ret*10+s[0]-'0';
    s.erase(0,1);
  }
  return ret;
}

int string_get_until_number(string &s)
{
  string_regularize(s, true);
  return string_get_number(s);
}

pair<pair<int,int>, int> string_get_node(string &s)
{
  string_regularize(s);
  while(s.length()>0 && s[0]!='n')
	s.erase(0,1);
  if (s.length()==0 || s[0]!='n')
    return make_pair(make_pair(-1,-1),-1);
  s.erase(0,1);
  int net=string_get_until_number(s);
  net_set.insert(net);
  int nodeX=string_get_until_number(s);
  int nodeY=string_get_until_number(s);
  auto ret=make_pair(make_pair(nodeX,nodeY), net);
  auto insert_ret=node_set.insert(ret);
  if (insert_ret.second){
    node_id[ret]=++node_cnt;
    father[node_id[ret]]=node_id[ret];
  }
  return ret;
}

int get_father(int id){
  if (father[id]==id)
    return id;
  int ret=get_father(father[id]);
  father[id]=ret;
  return ret;
}

int main(int argc, char ** argv)
// input #1: spice file name
// input #2: solution file name
// input #3: data output name
{
  if (argc <= 3)
    return 0;

  string grid_filename(argv[1]);
  fstream grid_stream(grid_filename.c_str(), ios_base::in);

  map<int,map<int,double>> grid_resist;
  map<int,map<int,double>> node_vias;
  map<int, double> node_package_voltage;
  map<int, double> node_package_resist;
  map<int, double> current_load;
  map<int, double> node_voltage;

  string line_string;
  int line_cnt=0;
  while(getline(grid_stream, line_string)){
    line_cnt++;
    string_regularize(line_string);
    if (!line_string.size())
      continue;

    if (line_string[0]=='.')
      continue;

    //comment: start with *
    if (line_string[0]=='*'){
      line_string.erase(0,1);
      string_regularize(line_string);
      //cout<<line_string<<endl;
      if (line_string[0]=='c') // circuit *******
	continue;
      if (line_string[0]=='v'){ //vias 
	string_regularize(line_string, true);
	int vias_a=string_get_number(line_string);
	string_regularize(line_string, true);
	int vias_b=string_get_number(line_string);
	cout<<"vias: "<<vias_a<<" "<<vias_b<<endl;
	continue;
      }
      if (line_string[0]=='l'){//layer
	continue;
      }
      continue;
    }
    
    //current load: start with i
    if (line_string[0]=='i'){
      while(line_string[0]!=' ')
	line_string.erase(0,1);
      line_string.erase(0,1);
      bool rev=false;
      if (line_string[0]=='0'){
	//i2312 0 n1_x_y rv
	rev=true;
      }else{
	//i2132 n1_x_y 0 rv
      }
      auto node=string_get_node(line_string);
      double rv;
      while(line_string.back()==' ' || line_string.back()=='\n' || line_string.back()=='\r') line_string.erase(line_string.size()-1);
      sscanf(line_string.substr(line_string.find_last_of(' ')).c_str(),"%lf",&rv);
      if (rev)
	rv=-rv;
      
      current_load[node_id[node]]+=rv;
      //cout<<node.first<<","<<node.second.first<<","<<node.second.second<<"="<<rv<<endl;
      continue;
    }

    //vias: start with V
    if (line_string[0]=='V'){
      while(line_string[0]!=' ')
	line_string.erase(0,1);
      auto node1=string_get_node(line_string);
      auto node2=string_get_node(line_string);
      double rv;
      sscanf(line_string.c_str(),"%lf",&rv);
      node_vias[node_id[node1]][node_id[node2]]=rv;
      node_vias[node_id[node2]][node_id[node1]]=rv;
      if (rv<eps)
	father[get_father(node_id[node1])]=node_id[node2];
      continue;
    }

    //Resistors: start with R
    if (line_string[0] == 'R'){
      line_string.erase(0,1);
      string_get_number(line_string);

      auto node1=string_get_node(line_string);
      auto node2=string_get_node(line_string);

      //      cout<<node1.first<<":"<<node1.second.first<<":"<<node1.second.second<<endl;
      //      cout<<node2.first<<":"<<node2.second.first<<":"<<node2.second.second<<endl;
      //      return 0;

      double rv;
      sscanf(line_string.c_str(),"%lf",&rv);
      
      grid_resist[node_id[node1]][node_id[node2]]=rv;
      grid_resist[node_id[node2]][node_id[node1]]=rv;

      if (rv<eps)
	father[get_father(node_id[node1])]=node_id[node2];

      continue;
    }

    //resistors of packages: start with r
    if (line_string[0] == 'r'){
      while(line_string[0]!=' ') 
	line_string.erase(0,1);
      auto node=string_get_node(line_string);
      double rv;
      while(line_string.back()==' ' || line_string.back()=='\n' || line_string.back()=='\r') line_string.erase(line_string.size()-1);
      sscanf(line_string.substr(line_string.find_last_of(' ')).c_str(),"%lf", &rv);
      //cout<<node.first<<":"<<node.second.first<<":"<<node.second.second<<"="<<rv<<endl;
      node_package_resist[node_id[node]]=rv;
      continue;
    }
    
    //voltage of packages: start with v
    if (line_string[0] == 'v'){
      while(line_string[0]!=' ')
	line_string.erase(0,1);
      auto node=string_get_node(line_string);
      double rv;
      while(line_string.back()==' ' || line_string.back()=='\n' || line_string.back()=='\r') line_string.erase(line_string.size()-1);
      sscanf(line_string.substr(line_string.find_last_of(' ')).c_str(),"%lf", &rv);
      node_package_voltage[node_id[node]]=rv;
      //cout<< node.first <<","<<node.second.first<<","<<node.second.second<<": "<<rv  <<endl;
      continue;
    }

  }

  cout<<"line cnt = "<<line_cnt<<endl;
  
  cout<<"current node cnt = "<<node_cnt<<endl;

  
  //-------------read solution---------------
  string solution_filename(argv[2]);
  fstream solution_stream(solution_filename.c_str(), ios_base::in);
  while(getline(solution_stream, line_string)){
    if (line_string[0]=='_'){
      // X_n3_xxxxx_xxxxx xxxx.xxxx
      continue;
    }
    auto node = string_get_node(line_string);
    if (node.first.first == -1){
      //maybe:
      // 1) G 0.000e+00
      continue;
    }
    double rv;
    sscanf(line_string.c_str(), "%lf", &rv);
    node_voltage[node_id[node]] = rv;
  }

  //----------set final ID------------------
  int final_node_cnt=0;
  for(int i=1; i<=node_cnt; i++){
    int ret=get_father(i);
    if (ret==i){
      final_id[i]=++final_node_cnt;
    }
  }
  for(int i=1; i<=node_cnt; i++){
    final_id[i]=final_id[get_father(i)];
  }
  cout<<"final node cnt = "<<final_node_cnt<<endl;
  
  // ------------write result-----------------

  string data_filename(argv[3]);
  fstream data_stream(data_filename.c_str(), ios_base::out);
  data_stream.precision(10);
  
  data_stream << net_set.size() << endl << endl;
  data_stream << final_node_cnt << endl << endl;

  map<int,map<int,double> > final_row;
  map<int, double> final_b;
  for(auto it: node_id){
    int id = it.second;

    //grid resistors
    for(auto it2: grid_resist[id]){
      int id2=it2.first;
      double rv=it2.second;
      if (rv<eps)
	continue;
      final_row[final_id[id]][final_id[id2]]+=-1.0f/rv;
      final_row[final_id[id]][final_id[id]]+=1.0f/rv;
    }

    //vias between nets
    for(auto it2: node_vias[id]){
      int id2=it2.first;
      double rv=it2.second;
      if (rv<eps)
	continue;
      cout<<"unexpected: vias not zero!"<<endl;
      final_row[final_id[id]][final_id[id2]]+=-1.0f/rv;
      final_row[final_id[id]][final_id[id]]+=1.0f/rv;
    }

    //current loads
    final_b[final_id[id]]+=-current_load[id];

    //package connections
    if (node_package_resist.find(id)!=node_package_resist.end()){
      final_row[final_id[id]][final_id[id]]+=1.0f/node_package_resist[id];
      final_b[final_id[id]]+=node_package_voltage[id]/node_package_resist[id];
    }
  }

  for(auto it: node_id){
    int id = it.second;
    if (get_father(id)!=id)
      continue;
    data_stream << final_id[id] << endl << it.first.second << " " << it.first.first.first << " " << it.first.first.second << endl;
    //output matrix row
    data_stream<<'\t'<<final_row[final_id[id]].size()<<endl;
    data_stream<<'\t';
    for (auto it2: final_row[final_id[id]])
      data_stream<<it2.first<<" "<<it2.second<<"  ";
    data_stream<<endl;
    data_stream<<'\t'<<final_b[final_id[id]]<<endl;
    data_stream<<node_voltage[id]<<endl;
    data_stream<<endl;
  }

  return 0;
}
