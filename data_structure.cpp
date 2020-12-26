#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <bits/stdc++.h>
#include <fstream>
#include <ctime>
#include <time.h>
#include <cmath>
#include "data_structure.h"

using namespace std;

/*
  Class Graph's methods' definition
*/
// Constructor
Graph::Graph() {

    nb_vertices = nb_edges = 0;
}

// Function used to add a new edge to the graph
void Graph::add_edge(struct edge_struct& newedge) {

    edges_id_tab.push_back(newedge.e_id);
    edges_internal_map[newedge.e_id] = newedge;

    ++nb_edges;

    // Get the vertices
    if(all_vertices.empty())
    {
        all_vertices.push_back(newedge.source);
        all_vertices.push_back(newedge.destination);
        //nb_vertices += 2;
    }
    else
    {
        if(!find_element_in_tab(all_vertices, newedge.source))
        {
            all_vertices.push_back(newedge.source);
            //++nb_vertices;
        }
        if(!find_element_in_tab(all_vertices, newedge.destination))
        {
            all_vertices.push_back(newedge.destination);
            //++nb_vertices;
        }
    }
    nb_vertices = all_vertices.size();

}

// Function used to add a vertices_map
void Graph::add_vertices_map(map<int, struct edge_endpoint>& vertices_internal_map_) {
    vertices_internal_map = vertices_internal_map_;
}

// Useful function to update an edge's residual bandwidth after a vlink treatment
void Graph::update_edge_bw_by_id(int const & i, long const & bandwidth_) {
    edges_internal_map[i].residual_bandwidth -= bandwidth_;
}

// Useful function to update a vertex's available cpu after a vnf placement
void Graph::update_vertex_cpu_by_id(int const & i, long const & cpu_) {
    vertices_internal_map[i].available_cpu = cpu_;
}

// We use this function to update a vertex's infos
void Graph::update_vertex_infos(int const & i, int const & cpu, int const & ram, int const & capacity) {
    vertices_internal_map[i].available_cpu -= cpu;
    vertices_internal_map[i].available_ram -= ram;
    vertices_internal_map[i].available_disk_capacity -= capacity;
}

// We just keep track of the really used nodes in the entire graph
void Graph::update_used_nodes_map(int const & i) {
    bool found(false);
    if(used_nodes_map.find(i) != used_nodes_map.end())
        found = true;

    if(!found)
        used_nodes_map[i] = true;
}

// We keep track of the really used edges in the graph
void Graph::update_used_edges_map(int const & i) {
    bool found(false);
    if(used_edges_map.find(i) != used_edges_map.end())
        found = true;

    if(!found)
        used_edges_map[i] = true;
}

// This function returns the number of vertices
int Graph::get_nb_vertices() const {
    return nb_vertices;
}

// This function returns the number of edges
int Graph::get_nb_edges() const {
    return nb_edges;
}

// This function returns a pointer on a specific subset
struct subset Graph::get_subset(int const & i) {
    return subsets[i];
};

// A useful recursive function to find root of a set
int Graph::find_root(struct subset subsets[], int i) {
	// find root and make root as parent of i (path compression)
	if (subsets[i].parent != i)
    {
		subsets[i].parent = find_root(subsets, subsets[i].parent);
    }

	return subsets[i].parent;
}

// The Union part of the Union-find algorithm is implemented right here
void Graph::union_subsets(struct subset subsets[], int x, int y) {
	int xroot(Graph::find_root(subsets, x));
	int yroot(Graph::find_root(subsets, y));

	// Attach smaller rank tree under root of high rank tree
	// (Union by Rank)
	if (subsets[xroot].rank < subsets[yroot].rank)
    {
		subsets[xroot].parent = yroot;
    }

	else if (subsets[xroot].rank > subsets[yroot].rank)
    {
 		subsets[yroot].parent = xroot;
    }


	// If ranks are same, then make one as root and increment
	// its rank by one
	else
	{
		subsets[yroot].parent = xroot;
		subsets[xroot].rank++;
	}
}

// Implementation of the Union-find process
bool Graph::union_find_traversal() {
    int V(nb_vertices), E(nb_edges);
    bool has_cycle(false);

	// Allocate memory for creating V sets
	subsets = (struct subset*) malloc( V * sizeof(struct subset) );
	//subsets = new subset();
	//struct subset subsets[V];
	for (int v(0); v < V; ++v)
	{
		subsets[v].parent = v;
		subsets[v].rank = 0;
	}

	// Iterate through all edges of graph, find sets of both
	// vertices of every edge, if sets are same, then there is
	// cycle in graph.
	for(int e(1); e <= E; ++e)
	{
		int x = Graph::find_root(subsets, edges_internal_map[e].source);
		int y = Graph::find_root(subsets, edges_internal_map[e].destination);

		if (x == y) has_cycle = true;

		Graph::union_subsets(subsets, x, y);
	}

    return has_cycle;

}

// We're getting all the subsets informations
string Graph::get_all_subsets_infos() {
    string result("All subsets informations :\n");
    int V(nb_vertices);
	for (int v(0); v < V; ++v)
	{
	    result += "\t====>Vertex " + to_string(v) + " ['parent':" + to_string(subsets[v].parent)
             + ", 'rank':" + to_string(subsets[v].rank)/* + ", 'path_to_root':" + to_string(get_vertices_path(v)) */
             + "]\n";

        vertex_with_parent_map[subsets[v].parent].push_back(v);
	}
    return result;

}

// An overview of a Graph object
void Graph::get_graph_recap() {
    // Declaration & Initialization
    string output("\nGraph recap /// Nb Edges => ");

    // Processing
    output += to_string(Graph::get_nb_edges()) + ", Nb Vertices => " + to_string(Graph::get_nb_vertices()) + ", ";
    (Graph::union_find_traversal()) ? output += "Graph contains a cycle" : output += "Graph DOES NOT contain a cycle";
    output += "\n" + Graph::get_all_subsets_infos();

    // Output
    cout << output << endl;
}

// We're getting the corresponding path, starting from a specific node inside the graph
map<int, vector<int>> Graph::get_path_from_vertex_map() const {
    return path_from_vertex_map;
}

// We can retrieve a specific node's struct, using its ID
struct edge_endpoint Graph::get_vertex_infos_by_id(int const & i) {
    return vertices_internal_map[i];
}

// We can retrieve a specific node's address, using its ID
struct edge_endpoint* Graph::get_vertex_addr_by_id(int const & i) {
    return &vertices_internal_map[i];
}

// We can retrieve a specific edge's address, using its ID
struct edge_struct* Graph::get_edge_addr_by_id(int const & i) {
    return &edges_internal_map[i];
}

// Getting all the edges IDs
vector<int> Graph::get_edges_id_tab() const {
    return edges_id_tab;
}

// We count the number of really used nodes inside the entire graph
unsigned int Graph::get_used_nodes_map_size() const {
    return used_nodes_map.size();
}

// Getting the used_nodes_map
map<int, bool> Graph::get_used_nodes_map() const {
    return used_nodes_map;
}

// We count the number of really used edges inside the entire graph
unsigned int Graph::get_used_edges_map_size() const {
    return used_edges_map.size();
}
// Getting the used_edges_map
map<int, bool> Graph::get_used_edges_map() const {
    return used_edges_map;
}

// This one is used just after the union_find_traversal process,
// so that we get all the possible paths to explore the entire graph
void Graph::check_all_possible_paths() {

    int V(nb_vertices);
    int parent_vertex;

    cout << "Checking all the paths starting from each vertex..." << endl;

	for (int v(0); v < V; ++v)
	{
        parent_vertex = subsets[v].parent;
        path_from_vertex_map[v].push_back(v);

        vector<int> result(vertex_with_parent_map.find(parent_vertex)->second);

        if(erase_element_in_tab(result, v))
        {
            if(v==0)
                sort(result.begin(), result.end());
            else
                reverse(result.begin(), result.end());

            for (auto x : result)
                path_from_vertex_map[v].push_back(x);
        }
        else
        {
            cout << "Execution exception!";
            break;
        }
	}

    cout << "\t====>Checking completed!" << endl;
}



/*
  Class VNC methods' definition
*/
// Constructor
VNC::VNC() {
    nb_vnf = nb_vlinks = 0;
}

// A function that returns the number of VNF
int VNC::get_nb_vnf() const {
    return nb_vnf;
}

// A function that returns the number of virtual links
int VNC::get_nb_vlinks() const {
    return nb_vlinks;
}

// This function is called to add a new VNF to the chain
void VNC::add_vnf(struct vnf& vnf_) {
    if (vnf_internal_map.find(vnf_.vnf_id)!= vnf_internal_map.end())
        cout << "VNF ID already exists! Please define a new one!\n";
    else
    {
        ++nb_vnf;
        vnf_internal_map[vnf_.vnf_id] = vnf_;
        vnf_id_tab.push_back(vnf_.vnf_id);
        vnf_placement_status[vnf_.vnf_id] = false;
    }
}

// This function helps with a virtual link adding between 2 vnf
void VNC::add_vlink(struct v_links& vlink_) {
    if (check_vlink_existence(vlink_.source_vnf, vlink_.destination_vnf))
        cout << "Virtual link already added! Please define a new one!\n";
    else
    {
        ++nb_vlinks;
        vlink_internal_map[vlink_.v_link_id] = vlink_;
        vlink_id_tab.push_back(vlink_.v_link_id);
        vlink_process_status.insert( pair<int, struct vlink_support>(vlink_.v_link_id, {vlink_.v_link_id,false,"","","",""}) );
    }
}

// We use this function to gather all the required informations after a vnf placement
void VNC::update_placement_result_map(struct vnf* pvnf, struct v_links* pvlink, struct edge_endpoint* pnode) {
    //placement_result_map.insert({pvnf, {pvlink, pnode}});
    placement_result_map[pvnf].insert({pvlink, pnode});
}

// We're saving the VNC infos' recap, for future purposes
void VNC::add_vnc_infos_map(int const & key, std::string const & given_value) {
    switch(key)
    {
    case 0:
        vnc_infos_map["vnf_infos"] = given_value;
        break;
    case 1:
        vnc_infos_map["vlink_infos"] = given_value;
        break;
    default:
        cout << "Error! Unknown key!" << endl;
    }
}

// We're updating the vnf placement status
void VNC::update_vnf_placement_status(int const & i, bool const & status) {
    if(vnf_placement_status.find(i) != vnf_placement_status.end())
        vnf_placement_status[i] = status;
    else
        cout << "Error! Unknown key!" << endl;
}

// We'll use this function to reset a VNC object t its initial values
void VNC::reset() {
    nb_vnf = nb_vlinks = 0;
    vnf_internal_map.clear();
    vlink_internal_map.clear();
    vnc_infos_map.clear();
    vnf_placement_status.clear();
    placement_result_map.clear();
    vnf_id_tab.clear();
    vlink_id_tab.clear();
}

// Getting all the VNF IDs
vector<int> VNC::get_vnf_id_tab() const {
    return vnf_id_tab;
}

// Getting all the virtual links IDs
vector<int> VNC::get_vlink_id_tab() const {
    return vlink_id_tab;
}

// Getting a VNF infos (address), using its ID
struct vnf* VNC::get_vnf_by_id(int const & i) {
    return &vnf_internal_map[i];
}

// Getting a virtual link infos (address), using its ID
struct v_links* VNC::get_vlink_by_id(int const & i) {
    return &vlink_internal_map[i];
}

// Getting the VNC recap infos inside vnc_recap_map
map<string, string> VNC::get_vnc_infos_map() const {
    return vnc_infos_map;
}

// We're getting the vnf_placement_status map and a specific value as well
map<int, bool> VNC::get_vnf_placement_status() const {
    return vnf_placement_status;
}
bool VNC::get_singlevnf_placement_status(int const & i) {
    return vnf_placement_status[i];
}

// We retrieve the vlink_process_status map and a specific value's address as well
map<int, struct vlink_support> VNC::get_vlink_process_status() const {
    return vlink_process_status;
}
struct vlink_support* VNC::get_vlink_support(int const & i) {
    return &vlink_process_status[i];
}

// Getting the VNC recap
string VNC::get_vnc_recap() {
    string recap_vnf("All VNF informations\n"), recap_vlink("");

    for(int i(0); i< nb_vnf; ++i)
    {
        recap_vnf += "\t====>VNF " + to_string(vnf_internal_map[i].vnf_id) + " ['description':" + vnf_internal_map[i].description
            + ", 'cpu':" + to_string(vnf_internal_map[i].cpu) + ", 'ram':" + to_string(vnf_internal_map[i].ram)
            + ", 'disk_capacity':" + to_string(vnf_internal_map[i].disk_capacity) + "]\n";
    }

    VNC::add_vnc_infos_map(0, recap_vnf);

    recap_vlink += "All VLinks informations\n";

    for(int i(0); i<nb_vlinks; ++i)
    {
        recap_vlink += "\t====>VLink " + to_string(vlink_internal_map[i].v_link_id) + " ['desc':" + vlink_internal_map[i].description
            + ", 'req_bw':" + to_string(vlink_internal_map[i].requested_bandwidth)
            + ", 'req_latency':" + to_string(vlink_internal_map[i].requested_latency)
            + ", 'vnf_src_id':" + to_string((vlink_internal_map[i].source_vnf)->vnf_id)
            + ", 'vnf_dest_id':" + to_string((vlink_internal_map[i].destination_vnf)->vnf_id) + "]\n";
    }

    VNC::add_vnc_infos_map(1, recap_vlink);

    return (recap_vnf+recap_vlink);
}

// With this function, we always avoid duplicated entries inside the vlink internal map
bool VNC::check_vlink_existence(struct vnf* source, struct vnf* dest) {
    bool already_exists(false);
    for (auto& it : vlink_internal_map) {
        if (it.second.source_vnf == source && it.second.destination_vnf == dest)
        {
            already_exists = true;
            break;
        }
    }
    return already_exists;
}



/*
  All the other methods' definition
*/
// Useful function to check the presence of an element i inside a tab
bool find_element_in_tab(vector<int> & tab, int const & i) {
    bool result(false);
    for(vector<int>::iterator it=tab.begin(); it!=tab.end(); ++it)
    {
        if(*it==i)
        {
            result = true;
            break;
        }
    }
    return result;
}

// Useful function to erase an element i inside a tab
bool erase_element_in_tab(vector<int> & tab, int const & i) {
    bool result(false);
    for(vector<int>::iterator it=tab.begin(); it!=tab.end(); ++it)
    {
        if(*it==i)
        {
            tab.erase(it);
            result = true;
            break;
        }
    }
    return result;
}

// Function used to initialize a edge_struct object
void edge_struct_init(struct edge_struct& edge_, int const & node_, long const & bandwidth_,
                        long const & residual_bandwidth_, int const & delay_, int const & cost_,
                        int const & source_, int const & destination_
) {
    edge_.e_id = node_;
    edge_.bandwidth = edge_.residual_bandwidth = bandwidth_;
    //edge_.residual_bandwidth = residual_bandwidth_;
    edge_.delay = delay_;
    edge_.cost = cost_;
    edge_.source = source_;
    edge_.destination = destination_;
}

// Function used to initialize a edge_endpoint object
void edge_endpoint_init(struct edge_endpoint& vertex_, int const & v_id_, int const & available_disk_capacity_,
                         int const & available_ram_, int const & available_cpu_, int const & cost_) {
    vertex_.v_id = v_id_;
    vertex_.available_disk_capacity = available_disk_capacity_;
    vertex_.available_ram = available_ram_;
    vertex_.available_cpu = available_cpu_;
    vertex_.cost = cost_;
}

// Checking path starting from a given vertex inside a graph
string check_path(Graph & mygraph, int const & i) {
    string path;
    vector<int> myvector(mygraph.get_path_from_vertex_map()[i]);

    for(auto x : myvector)
        path += " " + to_string(x);

    return path;
}

// Useful functions to sort all the vertices using the CPU and RAM values
// The vertex with the highest values will be choosed, and the next one...
bool comparator(pair<int, struct edge_endpoint>& a, pair<int, struct edge_endpoint>& b) {
    //return (a.second.available_cpu >= b.second.available_cpu && a.second.available_ram >= b.second.available_ram);
    return (a.second.available_cpu >= b.second.available_cpu && a.second.available_ram >= b.second.available_ram && a.second.available_disk_capacity >= b.second.available_disk_capacity);
}
vector<int> sort_vertices_by_resources(map<int, struct edge_endpoint>& vmap) {
    vector<int> sorted_vertices;

    // Vector of pairs' declaration
    vector<pair<int, struct edge_endpoint> > A;

    // Copy key-value pair from vmap to vector of pairs
    for (auto& it : vmap) {
        A.push_back(it);
    }

    // Sort using comparator function
    sort(A.begin(), A.end(), comparator);

    // Print the sorted value
    for (auto& it : A) {

        //cout << it.first << ' ' << it.second.bandwidth << endl;
        sorted_vertices.push_back(it.second.v_id);
    }

    return sorted_vertices;
}

// This function retrieves the vnf placement result and saves it into a file
void save_vnfplacement_result(string const & result_) {
    // Declarations & Initializations
    string const sys_datetime(get_sys_datetime());
    string const path_to_report("reports/VNFReport_"+sys_datetime +".txt");
    string const path_to_log("log/log_ind.txt");
    bool report_created(true);
    string mylog_output("[" + sys_datetime + "] ");

    // Opening the files stream
    ofstream mylog_stream(path_to_log.c_str(), ios::app); // We open the file and append new data to the existing
    ofstream myreport_stream(path_to_report.c_str());

    // Processing
    if(myreport_stream)
    {
        myreport_stream << result_ << endl;
    }
    else
    {
        report_created = false;
    }

    mylog_output += (report_created) ? "Report successfully created at " + path_to_report : "Report failure ! Something went wrong";
    mylog_stream << mylog_output << endl;

}

// We're getting the local system date-time and retrieve it as a string
string get_sys_datetime() {
    string sys_datetime(""),sys_year,sys_month,sys_day,sys_hour,sys_minute,sys_second;
    time_t timer1;
    time(&timer1);
    struct tm *newTime1;
    newTime1 = localtime(&timer1);
    sys_year = to_string(newTime1->tm_year+1900);
    sys_month = ((newTime1->tm_mon+1)<10) ? '0' + to_string(newTime1->tm_mon+1) : to_string(newTime1->tm_mon+1);
    sys_day = (newTime1->tm_mday<10) ? '0' + to_string(newTime1->tm_mday) : to_string(newTime1->tm_mday);
    sys_hour = (newTime1->tm_hour<10) ? '0' + to_string(newTime1->tm_hour) : to_string(newTime1->tm_hour);
    sys_minute = (newTime1->tm_min<10) ? '0' + to_string(newTime1->tm_min) : to_string(newTime1->tm_min);
    sys_second = (newTime1->tm_sec<10) ? '0' + to_string(newTime1->tm_sec) : to_string(newTime1->tm_sec);

    sys_datetime += sys_year+sys_month+sys_day+'_'+sys_hour+sys_minute+sys_second;

    return sys_datetime;
}

// We're using this function to format any string based on a field key's length
string format_string(int & field_length, string & field_value) {
    int fv_length(field_value.size()), mark;
    mark = field_length - fv_length;
    for(int i(0);i<mark;++i)
    {
        field_value += " ";
    }

    return field_value;
}

// Useful function to retrieve the index of an element in a vector
int get_elt_index(vector<int>& vec, int const & elt)
{
    auto it = find(vec.begin(), vec.end(), elt);
    int index(0);

    if (it != vec.end())
    {
        index = it - vec.begin();
    }
    else
    {
        index = -1;
    }
    return index;
}

// Useful function to find a convenient node for a specific vnf
int find_convenient_node_for_vnf(struct vnf* vnf_, Graph & graph_, vector<int> & path_) {
    int node_id;
    bool node_found(false);

    for(auto& node : path_)
    {
        if(vnf_->cpu <= graph_.get_vertex_infos_by_id(node).available_cpu && vnf_->ram <= graph_.get_vertex_infos_by_id(node).available_ram
            && vnf_->disk_capacity <= graph_.get_vertex_infos_by_id(node).available_disk_capacity)
        {
            graph_.update_vertex_infos(node, vnf_->cpu, vnf_->ram, vnf_->disk_capacity);
            node_found = true;
            node_id = node;
            break;
        }
    }

    //return (node_found) ? node_id : -1;
    if(node_found)
    {
        graph_.update_used_nodes_map(node_id);
        return node_id;
    }
    else
        return -1;
}

// We use this function to check if a virtual link's capacities can be supported by a specified edge in the graph
bool vlink_supported_by_edge(Graph & graph_, VNC* vnc_, int const & vlink_, int const & edge_src, int const & edge_dst) {
    bool result(false);
    long requested_bandwidth(vnc_->get_vlink_by_id(vlink_)->requested_bandwidth);
    vector<int> edges_id_tab = graph_.get_edges_id_tab();

    for(auto& id : edges_id_tab)
    {
        if((graph_.get_edge_addr_by_id(id)->source == edge_src && graph_.get_edge_addr_by_id(id)->destination == edge_dst)
           || (graph_.get_edge_addr_by_id(id)->source == edge_dst && graph_.get_edge_addr_by_id(id)->destination == edge_src))
        {
            //if(graph_.get_edge_addr_by_id(id)->bandwidth >= requested_bandwidth)
            if(graph_.get_edge_addr_by_id(id)->residual_bandwidth >= requested_bandwidth)
            {
                graph_.update_edge_bw_by_id(id, requested_bandwidth);
                graph_.update_used_edges_map(id);
                result = true;
                break;
            }
        }
    }

    return result;
}

// Final version of the VNF placement inside a given graph
// For that one, we check whether the graph's total bw and cpu can handle an entire VNC or not
// We verify each node's capacities as well, so that we can place the convenient vnf inside of it
bool vnf_placement_in_graph(Graph & mygraph_, vector<int> & sorted_vertices, VNC* myvnc, struct report_recap* recap_ptr) {
//string vnf_placement_in_graph(Graph & mygraph_, vector<int> & sorted_vertices, map<int, VNC> & vnc_map) {
    // Declarations & Initializations
    int nb_vertices(mygraph_.get_nb_vertices()), nb_edges(mygraph_.get_nb_edges());
    int nb_accepted_vnf(0), nb_rejected_vnf(0), nb_used_nodes(0), nb_used_edges(0);
    int real_nb_accepted_vnf(0), real_nb_rejected_vnf(0);
    int nb_vnf(myvnc->get_nb_vnf()), node_vnf_src, node_vnf_dst;
    int vldesc_ln(20), srcvnf_ln(29), dstvnf_ln(29), glk_ln(83);
    string final_output(""), output(""), final_path, output_if_succeeded, vlink_support_recap, output_if_failed("VNC placement failed!");
    string tmp_str1, nodes_recap(""), metrics_recap("Collected metrics\n"), nodes_resources_recap(""), edges_resources_recap(""), rejected_vnf_recap("");
    bool accepted(true), vnc_placed;
    vector<int> vnf_id_tab = myvnc->get_vnf_id_tab();
    struct vnf* vnf_src;
    struct vnf* vnf_dst;
    //struct v_links* tmp_vlink;
    map<int, int> vnf_placed_map;
    clock_t process_start = clock(), process_stop;

    // Major mapping process
    for(vector<int>::iterator it=sorted_vertices.begin(); it!= sorted_vertices.end(); ++it)
    {
        // Declarations & Initializations
        nb_accepted_vnf = nb_rejected_vnf = 0;
        Graph mygraph = mygraph_;
        final_path = "Accepted final path using vertices ID : ";
        output_if_succeeded = "[\n  ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------";
        output_if_succeeded += "\n  | VLink ID | VLink Description  | Src VNF --> Node in graph   | Dest VNF --> Node in graph  | Edges used in the graph                                                           |" ;
        output_if_succeeded += "\n  ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n";
        vlink_support_recap = output_if_succeeded;
        nodes_recap = "All nodes informations (Before placement)\n";
        accepted = true;
        vnf_placed_map.clear();

        // Getting the path starting from the vertex *it
        vector<int> mypath = mygraph.get_path_from_vertex_map()[*it];

       // Serializing the valid path
       for(vector<int>::iterator it0=mypath.begin(); it0!=mypath.end(); ++it0)
       {
            final_path += to_string(*it0);
            if ((it0+1)!=mypath.end())
                final_path += " -> ";

            nodes_recap += "\t====>Node " + to_string((mygraph.get_vertex_addr_by_id(*it0))->v_id)
                + " ['cpu':" + to_string((mygraph.get_vertex_addr_by_id(*it0))->available_cpu) + ", 'ram':"
                + to_string((mygraph.get_vertex_addr_by_id(*it0))->available_ram) + ", 'disk_capacity':"
                + to_string((mygraph.get_vertex_addr_by_id(*it0))->available_disk_capacity) + "]\n";
       }

        // We will iterate through the VNC virtual links
        for(auto& link_id : myvnc->get_vlink_id_tab())
        {
            vnf_src = myvnc->get_vlink_by_id(link_id)->source_vnf;
            vnf_dst = myvnc->get_vlink_by_id(link_id)->destination_vnf;
            struct vlink_support* vlink_ptr = myvnc->get_vlink_support(link_id);


            tmp_str1 = " " + myvnc->get_vlink_by_id(link_id)->description + " ";
            string mytmp_str2 = "  | VLink " + to_string(link_id) + "  |" + format_string(vldesc_ln, tmp_str1) + "|";
            vlink_ptr->description = mytmp_str2;
            output_if_succeeded += mytmp_str2;

            if (vnf_placed_map.find(vnf_src->vnf_id) == vnf_placed_map.end())
            {
                node_vnf_src = find_convenient_node_for_vnf(vnf_src, mygraph, mypath);
                if(node_vnf_src == -1)
                {
                    accepted = false;
                    ++nb_rejected_vnf;
                    break;
                }
                else
                {
                    vnf_placed_map[vnf_src->vnf_id] = node_vnf_src;
                    myvnc->update_placement_result_map(vnf_src, myvnc->get_vlink_by_id(link_id), mygraph.get_vertex_addr_by_id(node_vnf_src));
                    //myvnc->update_vnf_placement_status(vnf_src->vnf_id, true);
                    ++nb_accepted_vnf;
                }
            }
            else
            {
                node_vnf_src = vnf_placed_map[vnf_src->vnf_id];
            }
            //tmp_str1 = " VNF " + to_string(vnf_src->vnf_id) + " ['desc':" + vnf_src->description + ", 'Placed in node':" + to_string(node_vnf_src) + "] ";
            tmp_str1 = " VNF " + to_string(vnf_src->vnf_id) + " --> Node " + to_string(node_vnf_src) + " ";
            mytmp_str2 = format_string(srcvnf_ln, tmp_str1) + "|";
            vlink_ptr->graph_src_node = mytmp_str2;
            output_if_succeeded += mytmp_str2;


            if (vnf_placed_map.find(vnf_dst->vnf_id) == vnf_placed_map.end())
            {
                node_vnf_dst = find_convenient_node_for_vnf(vnf_dst, mygraph, mypath);
                if(node_vnf_dst == -1)
                {
                    accepted = false;
                    ++nb_rejected_vnf;
                    break;
                }
                else
                {
                    vnf_placed_map[vnf_dst->vnf_id] = node_vnf_dst;
                    myvnc->update_placement_result_map(vnf_dst, myvnc->get_vlink_by_id(link_id), mygraph.get_vertex_addr_by_id(node_vnf_dst));
                    //myvnc->update_vnf_placement_status(vnf_dst->vnf_id, true);
                    ++nb_accepted_vnf;
                }
            }
            else
            {
                node_vnf_dst = vnf_placed_map[vnf_dst->vnf_id];
            }
            //tmp_str1 = " VNF " + to_string(vnf_dst->vnf_id) + " ['desc':" + vnf_dst->description + ", 'Placed in node':" + to_string(node_vnf_dst) + "] ";
            tmp_str1 = " VNF " + to_string(vnf_dst->vnf_id) + " --> Node " + to_string(node_vnf_dst) + " ";
            mytmp_str2 = format_string(dstvnf_ln, tmp_str1) + "|";
            vlink_ptr->graph_dest_node = mytmp_str2;
            output_if_succeeded += mytmp_str2;

            // We're then checking the graph's edge capacities
            if(node_vnf_src == node_vnf_dst)
            {
                tmp_str1 = " N/A: Same node ";
                vlink_ptr->supported = true;
                myvnc->update_vnf_placement_status(vnf_src->vnf_id, true);
                myvnc->update_vnf_placement_status(vnf_dst->vnf_id, true);
                //output_if_succeeded += format_string(glk_ln, tmp_str1) + "|\n";
            }
            else
            {
                //tmp_str1 = "(" + to_string(node_vnf_src) + "," + to_string(node_vnf_dst) + ") ";

                //tmp_str1 = " " + to_string(node_vnf_src) + '-' + to_string(node_vnf_dst) + " ";
                //output_if_succeeded += format_string(glk_ln, tmp_str1) + "|\n";

                map<int, vector<int>> vlink_inside_graph;
                tmp_str1 = " ";
                int real_node_src(0), real_node_dst(0);

                // Retrieving the nodes's indexes and taking actions
                if(get_elt_index(mypath,node_vnf_src) < get_elt_index(mypath,node_vnf_dst))
                {
                    real_node_src = node_vnf_src;
                    real_node_dst = node_vnf_dst;
                }
                else
                {
                    real_node_src = node_vnf_dst;
                    real_node_dst = node_vnf_src;
                }

                for(unsigned int i(0);i<mypath.size();++i)
                {
                    if(mypath[i] == real_node_src)
                    {
                        int ct(0);
                        unsigned int j(i);
                        do
                        {
                            vlink_inside_graph[ct].push_back(mypath[j]);
                            vlink_inside_graph[ct].push_back(mypath[j+1]);
                            //tmp_str1 += to_string(mypath[j]) + "-";
                            ++j;
                            if(mypath[j]==real_node_dst)
                                break;

                            ++ct;
                        } while(mypath[j]!=real_node_dst);
                    }

                }
                //tmp_str1 += to_string(node_vnf_dst);

                for(map<int,vector<int>>::iterator it1=vlink_inside_graph.begin();it1!=vlink_inside_graph.end();++it1)
                {
                    //bool tmp_link_supported(true);
                    if(vlink_supported_by_edge(mygraph, myvnc, link_id, it1->second[0], it1->second[1]))
                    {
                        tmp_str1 += "(" + to_string(it1->second[0]) + "," + to_string(it1->second[1]) + ") ";
                        vlink_ptr->supported = true;
                    }
                    else
                    {
                        nb_accepted_vnf -= 1;
                        //tmp_link_supported = false;
                        //if(myvnc->get_singlevnf_placement_status(vnf_src->vnf_id))
                        //    myvnc->update_vnf_placement_status(vnf_src->vnf_id, false);

                        //if(myvnc->get_singlevnf_placement_status(vnf_dst->vnf_id))
                        //    myvnc->update_vnf_placement_status(vnf_dst->vnf_id, false);
                        //accepted = false;
                        //output_if_failed = "Failure : At least one virtual link's bandwidth is not supported by a graph's edge!";
                        tmp_str1 += "";
                        break;
                    }

                }


            }
            mytmp_str2 = format_string(glk_ln, tmp_str1) + "|\n";
            vlink_ptr->graph_edges = mytmp_str2;
            output_if_succeeded += mytmp_str2;


        }


        if(accepted)
        {
            process_stop = clock();

			nodes_recap += "All nodes informations (After placement)\n";
			for(vector<int>::iterator it2=mypath.begin(); it2!=mypath.end(); ++it2)
			{
				nodes_recap += "\t====>Node " + to_string((mygraph.get_vertex_addr_by_id(*it2))->v_id)
					+ " ['cpu':" + to_string((mygraph.get_vertex_addr_by_id(*it2))->available_cpu) + ", 'ram':"
					+ to_string((mygraph.get_vertex_addr_by_id(*it2))->available_ram) + ", 'disk_capacity':"
					+ to_string((mygraph.get_vertex_addr_by_id(*it2))->available_disk_capacity) + "]\n";
			}


            for(auto& vl : myvnc->get_vlink_process_status())
            {
                if(vl.second.supported)
                {
                    vlink_support_recap += vl.second.description + vl.second.graph_src_node
                        + vl.second.graph_dest_node + vl.second.graph_edges;

                    vnf_src = myvnc->get_vlink_by_id(vl.second.v_id)->source_vnf;
                    vnf_dst = myvnc->get_vlink_by_id(vl.second.v_id)->destination_vnf;

                    if(!myvnc->get_singlevnf_placement_status(vnf_src->vnf_id))
                    {
                        myvnc->update_vnf_placement_status(vnf_src->vnf_id, true);
                    }
                    if(!myvnc->get_singlevnf_placement_status(vnf_dst->vnf_id))
                    {
                        myvnc->update_vnf_placement_status(vnf_dst->vnf_id, true);
                    }
                }
            }
            vlink_support_recap += "]";


            real_nb_accepted_vnf = real_nb_rejected_vnf = 0;
			for(auto& ct : myvnc->get_vnf_placement_status())
            {
                if(ct.second)
                    ++real_nb_accepted_vnf;
                else
                    ++real_nb_rejected_vnf;
            }



            if(real_nb_accepted_vnf==nb_vnf)
            {
                output = "VNC placement succeeded! Find below the recap\n";
            }
            else
            {
                output = "VNC placement partially complete! Find below the recap\n";

                map<int, bool> vnf_placement_status = myvnc->get_vnf_placement_status();
                //rejected_vnf_recap = "\t\t|\t\t Rejected VNF:\n";
                for(auto& v : vnf_placement_status)
                {
                    if(v.second == false)
                    {
                        rejected_vnf_recap += "\t\t|\t\t\t\t\t             ===> " + to_string(v.first) + " ['desc':"
                            + myvnc->get_vnf_by_id(v.first)->description + "]\n";
                    }

                }

            }

            final_path += '\n';
            output_if_succeeded += "]";


            nb_used_nodes = mygraph.get_used_nodes_map_size();
            nb_used_edges = mygraph.get_used_edges_map_size();

            for(auto& it2 : mygraph.get_used_nodes_map())
            {
                int residual_cpu((mygraph.get_vertex_addr_by_id(it2.first))->available_cpu);
                int residual_ram((mygraph.get_vertex_addr_by_id(it2.first))->available_ram);
                int residual_disk((mygraph.get_vertex_addr_by_id(it2.first))->available_disk_capacity);
                int initial_cpu((mygraph_.get_vertex_infos_by_id(it2.first)).available_cpu);
                int initial_ram((mygraph_.get_vertex_infos_by_id(it2.first)).available_ram);
                int initial_disk((mygraph_.get_vertex_infos_by_id(it2.first)).available_disk_capacity);
                int rate_used_cpu(0), rate_used_ram(0), rate_used_disk(0);

                rate_used_cpu = ceil(((initial_cpu-residual_cpu) * 100) / initial_cpu);
                rate_used_ram = ceil(((initial_ram-residual_ram) * 100) / initial_ram);
                rate_used_disk = ceil(((initial_disk-residual_disk) * 100) / initial_disk);


                nodes_resources_recap += "  \t\t|\t\t    ==>%Used resources Server " + to_string(it2.first)
                    + " ['cpu_used_at':" + to_string(rate_used_cpu) + "%, 'ram_used_at':" + to_string(rate_used_ram)
                    + "%, 'disk_capacity_used_at':" + to_string(rate_used_disk) + "%]\n";

            }

            for(auto& it3 : mygraph.get_used_edges_map())
            {
                long residual_bw((mygraph.get_edge_addr_by_id(it3.first))->residual_bandwidth);
                long initial_bw((mygraph.get_edge_addr_by_id(it3.first))->bandwidth);
                int rate_used_bw(0);

                rate_used_bw = ceil(((initial_bw-residual_bw) * 100) / initial_bw);

                edges_resources_recap += "  \t\t|\t\t       ==>%Used resources Edge " + to_string(it3.first)
                    + " ['bandwidth_used_at':" + to_string(rate_used_bw) + "%]\n";
            }

            mygraph_ = mygraph;

            break;
        }

    }
    //clock_t process_stop = clock();


    // Saving & Returning result
	//final_output += (accepted) ? (output+final_path+output_if_succeeded) : output_if_failed;
	if(accepted)
    {
        // Collecting the metrics
        double elapsed_time(0);
        int accepted_vnf_rate(0), rejected_vnf_rate(0), rate_used_nodes(0), rate_used_edges(0);

        accepted_vnf_rate = real_nb_accepted_vnf * 100 / nb_vnf;
        rejected_vnf_rate = 100 - accepted_vnf_rate;
        //rejected_vnf_rate = real_nb_rejected_vnf * 100 / nb_vnf;
        rate_used_nodes = ceil(nb_used_nodes * 100 / nb_vertices);
        rate_used_edges = ceil(nb_used_edges * 100 / nb_edges);

        elapsed_time = ((double)(process_stop-process_start)) / CLOCKS_PER_SEC;

        // Recap's construction
        metrics_recap += "  \t\t----------------------------------------------------------------------\n  \t\t| Accepted VNF (%)          \t\t: "
            + to_string(accepted_vnf_rate) + "\n  \t\t|\n"
            + "  \t\t| Rejected VNF (%)          \t\t: " + to_string(rejected_vnf_rate) + "  \t\t\n" + rejected_vnf_recap + "  \t\t|\n"
            + "  \t\t| Process elapsed CPU time (seconds): " + to_string(elapsed_time) + "\n  \t\t|\n"
            + "  \t\t| Nb really used servers    \t\t: " + to_string(nb_used_nodes) +  " / " + to_string(nb_vertices) + " servers ("
            + to_string(rate_used_nodes) +  "% used)\n" + nodes_resources_recap + "  \t\t|\n"
            + "  \t\t| Nb really used edges      \t\t: " + to_string(nb_used_edges) + " / " + to_string(nb_edges) + " edges ("
            + to_string(rate_used_edges) +  "% used)\n" + edges_resources_recap
            + "  \t\t----------------------------------------------------------------------\n";

        /*
        string report(final_output);
        report += "\t\t--------------------------------------------------------------\n\t\t    " + output
            + "  \t\t--------------------------------------------------------------\n"
            + "----------------------------------------------------------------------------------------\n"
            + metrics_recap + "----------------------------------------------------------------------------------------\n"
            + final_path + "----------------------------------------------------------------------------------------\n"
            + nodes_recap + "----------------------------------------------------------------------------------------\n"
            + myvnc->get_vnc_infos_map()["vnf_infos"] + "----------------------------------------------------------------------------------------\n"
            + myvnc->get_vnc_infos_map()["vlink_infos"] + "----------------------------------------------------------------------------------------\n\n"
            + output_if_succeeded;
        final_output += output+final_path+output_if_succeeded+"\n\n"+metrics_recap;

        // Creating the report file
        save_vnfplacement_result(report);
        */

        // Saving the recaps
        recap_ptr->processed = true;
        recap_ptr->elapsed_time = elapsed_time;
        recap_ptr->metrics_recap = metrics_recap;
        recap_ptr->output = output;
        recap_ptr->final_path = final_path;
        recap_ptr->nodes_recap = nodes_recap;
        recap_ptr->vnf_infos = myvnc->get_vnc_infos_map()["vnf_infos"];
        recap_ptr->vlinks_infos = myvnc->get_vnc_infos_map()["vlink_infos"];
        //recap_ptr->placement_result = output_if_succeeded;
        recap_ptr->placement_result = vlink_support_recap;

        vnc_placed = true;

    }
    else
    {
        //final_output += output_if_failed;
        vnc_placed = false;
    }


	//return final_output;
	return vnc_placed;
}

// This function will use the vnc_placement_in_graph() function to iterate through
// the entire vnc_map content
string requests_placement_process(Graph & mygraph_, vector<int> & sorted_vertices,
                                  std::map<int, VNC> & vnc_map, map<int, struct report_recap> & trequests_recap_map
) {
    string final_output(""),req_metrics_recap("Collected metrics\n"),req_recap_intro("Request(s) processed! Find below the recap\n");
    string nodes_resources_recap(""), edges_resources_recap(""), each_vnc_report_recap(""), tmp_vnc_recap(""), report("");
    int nb_nodes(mygraph_.get_nb_vertices()), nb_edges(mygraph_.get_nb_edges());
    int nb_accepted_req(0), nb_requests = (int) vnc_map.size(), nb_used_nodes(0), nb_used_edges(0);
    int accepted_req_rate(0), rejected_req_rate(0),rate_used_nodes(0), rate_used_edges(0);
    double proc_elapsed_time(0);
    Graph graph_init = mygraph_;
    vector<string> each_vnc_report;

    for(auto& vnc : vnc_map)
    {
        tmp_vnc_recap = vnc.second.get_vnc_recap();
        trequests_recap_map.insert( pair<int, struct report_recap>(vnc.first, {vnc.first,false,0,"","","","","","",""}) );
        //trequests_recap_map[vnc.first] = {vnc.first,false,"","","","","",""};
        struct report_recap* recap_ptr = &trequests_recap_map[vnc.first];
        VNC* vnc_ptr = & vnc.second;
        if(vnf_placement_in_graph(mygraph_, sorted_vertices, vnc_ptr, recap_ptr))
        {
            ++nb_accepted_req;
        }
    }

    // Requests placement metrics
    nb_used_nodes = mygraph_.get_used_nodes_map_size();
    nb_used_edges = mygraph_.get_used_edges_map_size();
    rate_used_nodes = ceil(nb_used_nodes *100 / nb_nodes);
    rate_used_edges = ceil(nb_used_edges * 100 / nb_edges);
    accepted_req_rate = nb_accepted_req * 100 / nb_requests;
    rejected_req_rate = 100 - accepted_req_rate;

    for(auto& re : trequests_recap_map)
    {
        proc_elapsed_time += re.second.elapsed_time;
        string output = "\t\t------------------------------------------------------------------------\n\t\t    Request #" + to_string(re.first);

        if(re.second.processed)
        {
            output += ": "+ re.second.output
                + "  \t\t------------------------------------------------------------------------\n"
                + "----------------------------------------------------------------------------------------\n"
                + re.second.metrics_recap + "----------------------------------------------------------------------------------------\n"
                + re.second.final_path + "----------------------------------------------------------------------------------------\n"
                //+ re.second.nodes_recap + "----------------------------------------------------------------------------------------\n"
                + re.second.vnf_infos + "----------------------------------------------------------------------------------------\n"
                + re.second.vlinks_infos + "----------------------------------------------------------------------------------------\n\n"
                + re.second.placement_result;
        }
        else
        {
            output += ": VNC placement failed!\n  \t\t--------------------------------------------------------------\n";
        }

        each_vnc_report.push_back(output);
    }

    for(auto& it : mygraph_.get_used_nodes_map())
    {
        int residual_cpu((mygraph_.get_vertex_addr_by_id(it.first))->available_cpu);
        int residual_ram((mygraph_.get_vertex_addr_by_id(it.first))->available_ram);
        int residual_disk((mygraph_.get_vertex_addr_by_id(it.first))->available_disk_capacity);
        int initial_cpu((graph_init.get_vertex_infos_by_id(it.first)).available_cpu);
        int initial_ram((graph_init.get_vertex_infos_by_id(it.first)).available_ram);
        int initial_disk((graph_init.get_vertex_infos_by_id(it.first)).available_disk_capacity);
        int rate_used_cpu(0), rate_used_ram(0), rate_used_disk(0);

        rate_used_cpu = ceil(((initial_cpu-residual_cpu) * 100) / initial_cpu);
        rate_used_ram = ceil(((initial_ram-residual_ram) * 100) / initial_ram);
        rate_used_disk = ceil(((initial_disk-residual_disk) * 100) / initial_disk);


        nodes_resources_recap += "  \t\t|\t        ==>%Used resources Server " + to_string(it.first)
            + " ['cpu_used_at':" + to_string(rate_used_cpu) + "%, 'ram_used_at':" + to_string(rate_used_ram)
            + "%, 'disk_capacity_used_at':" + to_string(rate_used_disk) + "%]\n";

    }

    for(auto& it1 : mygraph_.get_used_edges_map())
    {
        long residual_bw((mygraph_.get_edge_addr_by_id(it1.first))->residual_bandwidth);
        long initial_bw((mygraph_.get_edge_addr_by_id(it1.first))->bandwidth);
        int rate_used_bw(0);

        rate_used_bw = ceil(((initial_bw-residual_bw) * 100) / initial_bw);

        edges_resources_recap += "  \t\t|\t          ==>%Used resources Edge " + to_string(it1.first)
            + " ['bandwidth_used_at':" + to_string(rate_used_bw) + "%]\n";
    }

    req_metrics_recap  += "  \t\t----------------------------------------------------------------------\n  \t\t| Accepted requests (%)          \t: "
            + to_string(accepted_req_rate) + " (" + to_string(nb_accepted_req) + " / " + to_string(nb_requests) + " requests )\n  \t\t|\n"
            + "  \t\t| Rejected requests (%)          \t: " + to_string(rejected_req_rate) + "\n  \t\t|\n"
            //+ "  \t\t\n" + rejected_vnf_recap + "  \t\t|\n"
            + "  \t\t| Process elapsed CPU time (seconds)    : " + to_string(proc_elapsed_time) + "\n  \t\t|\n"
            + "  \t\t| Nb really used servers    \t\t: " + to_string(nb_used_nodes) +  " / " + to_string(nb_nodes) + " servers ("
            + to_string(rate_used_nodes) +  "% used)\n" //+ "  \t\t|\n"
            + nodes_resources_recap + "  \t\t|\n"
            + "  \t\t| Nb really used edges      \t\t: " + to_string(nb_used_edges) + " / " + to_string(nb_edges) + " edges ("
            + to_string(rate_used_edges) +  "% used)\n"
            + edges_resources_recap
            + "  \t\t----------------------------------------------------------------------\n";

    for(vector<string>::iterator rep=each_vnc_report.begin();rep!=each_vnc_report.end();++rep)
    {
        //string tmp_str = (*rep);
        each_vnc_report_recap += "\n" + (*rep) + "\n\n";
    }


    final_output += "\t\t--------------------------------------------------------------\n\t\t          " + req_recap_intro
        + "  \t\t--------------------------------------------------------------\n" + req_metrics_recap ;


    report = final_output + each_vnc_report_recap;

    // Creating the report file
    save_vnfplacement_result(report);

    return final_output;
}

// A useful function to load an entire infrastructure and the requests as well at once
void infra_and_requests_mass_loading(map<int, struct edge_struct>& edges_map, map<int, struct edge_endpoint>& vertices_map,
                                        //map<int, struct v_links>& vlinks_map, map<int, struct vnf>& vnf_map,
                                        Graph & mygraph, map<int, VNC>& vnc_map
) {
    /*
        We read the required informations from 2 given files and treat them to gather at once
        the nodes and the links inside a network, and the vnf and virtual links of a
        virtual network chain
    */
    cout << "\nCreation of the entire network infrastructure and the VNC elements..." << endl;

    // Infrastructure's treatment
    string const infra_path("input_files/infra.txt");
    string text_infra("");
    ifstream f_infra(infra_path.c_str());

    if(f_infra)
    {
        getline(f_infra, text_infra);
        //f_infra.ignore();

        vector<string> infra_vec = split_string_genesis(text_infra, 1);

        // Load of the network's nodes
        nodes_mass_loading(vertices_map, infra_vec[0]);
        mygraph.add_vertices_map(vertices_map);

        // Load of the network's physical links
        phys_links_mass_loading(edges_map, infra_vec[1]);
        for(auto& ed : edges_map)
        {
            mygraph.add_edge(ed.second);
        }

        //f_infra.close();
    }
    else
        cout << "\nError : File opening failure" << endl;


    // VNC requests' treatment
    string const requests_path("input_files/requests.txt");
    string text_requests("");
    ifstream f_req(requests_path.c_str());
    //VNC myvnc;

    if(f_req)
    {
        int iter(0);
        text_requests = "";
        while(getline(f_req, text_requests))
        {
            //myvnc.reset();
            //vnf_map.clear();
            //vlinks_map.clear();
            //VNC myvnc;
            vnc_map[iter] = VNC();
            map<int, struct v_links> vlinks_map;
            map<int, struct vnf> vnf_map;
            VNC* vnc_addr = &vnc_map[iter];

            vector<string> requests_vec = split_string_genesis(text_requests, 2);

            // Load of the virtual network functions
            vnf_mass_loading(iter, vnf_map, requests_vec[0]);
            for(auto& v : vnf_map)
            {
                //myvnc.add_vnf(v.second);
                vnc_addr->add_vnf(v.second);
            }

            // Load of the VNC's links
            vnc_links_mass_loading(iter, vnc_addr, vlinks_map, requests_vec[1]);
            for(auto& vl : vlinks_map)
            {
                //myvnc.add_vlink(vl.second);
                vnc_addr->add_vlink(vl.second);
            }

            //vnc_map.insert(pair<int, VNC>(iter, myvnc));
            //vnc_map[iter] = myvnc;

            ++iter;
            //f_req.ignore();
        }

        //f_req.close();

    }
    else
        cout << "\nError : File opening failure" << endl;

}

// The five next functions are used for the treatment of the input files's content
vector<string> split_string_genesis(const string &text, const int &key) {
    string sep1(""), sep2(""), text_part1(""), text_part2("");

    if(key==1)
    {
        // Infrastructure file
        sep1 = "NODES=";
        sep2 = "LINKS=";
    }
    else if(key==2)
    {
        // Requests file
        sep1 = "VNF=";
        sep2 = "VLINKS=";
    }
    else
    {
        sep1 = sep2 = "";
    }
    if(sep1=="" || sep2=="")
    {
        text_part1 = text_part2 = "N/A";
    }
    else
    {
        int ln_sep1 = (int) sep1.size(), ln_sep2 = (int) sep2.size(), ln_text = (int) text.size(), pos1 = text.find(sep2, 0);
        text_part1 = text.substr(ln_sep1, pos1-ln_sep1);
        text_part2 = text.substr(pos1+ln_sep2, ln_text);
    }

    return std::vector<string>{text_part1, text_part2};

}
void split_string(map<string, string> &mymap, const string &text, char sep) {
    int start = 0, end1 = 0;
    int count(0);
    while ((end1 = text.find(sep, start)) != (int) string::npos)
    {
        string key = to_string(count);
        string val = text.substr(start, end1 - start);
        mymap.insert(pair<string,string>(key, val));
        start = end1 + 1;
        ++count;
    }
}
void split_string_v2(map<string, string> &mymap, const string &text, char sep) {
    int start = 0, end1 = 0, end2 = (int) text.size();
    end1 = text.find(sep, start);
    string key = text.substr(start, end1);
    string val = text.substr(end1+1, end2-1);
    mymap.insert(pair<string,string>(key, val));
}
void split_string_links_v2(map<vector<string>, string> &mymap, const string &text, char sep) {
    int start = 0, end1 = 0, end2 = (int) text.size();
    end1 = text.find(sep, start);
    string key = text.substr(start, end1);
    string val = text.substr(end1+1, end2-1);

    int startv2 = 0, end1v2 = 0;
    //int end2v2 = (int) key.size();
    end1v2 = key.find(',', startv2);
    string key1 = key.substr(startv2+1, end1v2-1);
    string key2 = key.substr(end1v2+1, 1);

    bool keys_found(false);
    for(auto & it : mymap)
    {
        if((it.first[0]==key1 &&it.first[1]==key2) || (it.first[0]==key2 &&it.first[1]==key1))
        {
            keys_found = true;
            break;
        }
    }

    if(!keys_found)
    {
        mymap.insert(pair<vector<string>,string>({key1,key2},val));
    }
}
void split_string_v3(vector<string> & myvec, string &text, char sep) {
    int start = 0, end1 = 0, end2 = (int) text.size();
    if(end2 > 0)
    {
        end1 = text.find(sep, start);
        myvec.push_back(text.substr(start, end1));
        string mytext = text.substr(end1+1, end2-1);
        split_string_v3(myvec, mytext, sep);
    }
}

// A mass loading treatment for the nodes of the graph
void nodes_mass_loading(map<int, struct edge_endpoint>& vertices_map, string &text) {
    map<string, string> mymap_tmp;
    map<string, string> mymap_tmp1;
    map<string, vector<string>> mymap_tmp2;
    split_string(mymap_tmp, text, ';');

    for(auto& it : mymap_tmp)
    {
        split_string_v2(mymap_tmp1, it.second, ':');
        //cout << it.first << " ==> " << it.second << endl;
    }
    //cout << endl;
    for(auto& it : mymap_tmp1)
    {
        split_string_v3(mymap_tmp2[it.first], it.second, ',');
        //cout << it.first << " ==> " << it.second << endl;
    }
    //cout << endl;
    for(auto& it : mymap_tmp2)
    {
        int node_id(stoi(it.first)), node_disk(stoi(it.second[0])), node_ram(stoi(it.second[1]));
        int node_cpu(stoi(it.second[2])), node_cost(stoi(it.second[3]));

        vertices_map.insert({node_id, {node_id,node_disk,node_ram,node_cpu,node_cost}});
    }

    //cout << "\n[" << vertices_map[0].v_id << " " << vertices_map[0].available_disk_capacity << "]" << endl;
    cout << "\t====> " << vertices_map.size() << " Node(s) entries successfully loaded" << endl;

}

// A mass loading treatment for the graph's physical links
void phys_links_mass_loading(map<int, struct edge_struct>& edges_map, string &text) {
    map<string, string> mymap_tmp;
    map<vector<string>, string> mymap_tmp1;
    map<vector<string>, vector<string>> mymap_tmp2;
    split_string(mymap_tmp, text, ';');

    for(auto& it : mymap_tmp)
    {
        split_string_links_v2(mymap_tmp1, it.second, ':');
        //cout << it.first << " ==> " << it.second << endl;
    }
    //cout << endl;
    for(auto& it : mymap_tmp1)
    {
        split_string_v3(mymap_tmp2[it.first], it.second, ',');
        //cout << '(' << it.first[0] << ',' << it.first[1] << ')' << " ==> " << it.second << endl;
    }
    //cout << endl;
    for(auto& it : mymap_tmp2)
    {
        int link_id = (edges_map.size() > 0) ? (int) edges_map.size() : 0;
        int link_delay(stoi(it.second[1])), link_cost(stoi(it.second[2]));
        int link_src(stoi(it.first[0])), link_dest(stoi(it.first[1]));
        long link_bw(stol(it.second[0]));

        edges_map.insert({link_id, {link_id,link_bw,link_bw,link_delay,link_cost,link_src,link_dest}});
    }


    cout << "\t====> " << edges_map.size() << " Network link(s) entries successfully loaded" << endl;

}

// A mass loading treatment for the VNF
void vnf_mass_loading(int const & i, map<int, struct vnf>& vnf_map, string &text) {
    map<string, string> mymap_tmp;
    map<string, string> mymap_tmp1;
    map<string, vector<string>> mymap_tmp2;
    split_string(mymap_tmp, text, ';');

    for(auto& it : mymap_tmp)
    {
        split_string_v2(mymap_tmp1, it.second, ':');
        //cout << it.first << " ==> " << it.second << endl;
    }
    //cout << endl;
    for(auto& it : mymap_tmp1)
    {
        split_string_v3(mymap_tmp2[it.first], it.second, ',');
        //cout << it.first << " ==> " << it.second << endl;
    }
    //cout << endl;
    for(auto& it : mymap_tmp2)
    {
        string vnf_desc(it.second[0]);
        int vnf_id(stoi(it.first)), vnf_cpu(stoi(it.second[1])), vnf_ram(stoi(it.second[2])), vnf_disk(stoi(it.second[3]));

        vnf_map.insert({vnf_id, {vnf_id,vnf_desc,vnf_cpu,vnf_ram,vnf_disk}});
    }

    //cout << "\n[" << vertices_map[0].v_id << " " << vertices_map[0].available_disk_capacity << "]" << endl;
    cout << "\t====> [Request" << i << "] " << vnf_map.size() << " VNF entries successfully loaded" << endl;

}

// A mass loading treatment for the virtual links of the VNC
void vnc_links_mass_loading(int const & i, VNC* myvnc, map<int, struct v_links>& vlinks_map, string &text) {
    map<string, string> mymap_tmp;
    map<vector<string>, string> mymap_tmp1;
    map<vector<string>, vector<string>> mymap_tmp2;
    split_string(mymap_tmp, text, ';');

    for(auto& it : mymap_tmp)
    {
        split_string_links_v2(mymap_tmp1, it.second, ':');
        //cout << it.first << " ==> " << it.second << endl;
    }
    //cout << endl;
    for(auto& it : mymap_tmp1)
    {
        split_string_v3(mymap_tmp2[it.first], it.second, ',');
        //cout << '(' << it.first[0] << ',' << it.first[1] << ')' << " ==> " << it.second << endl;
    }
    //cout << endl;
    for(auto& it : mymap_tmp2)
    {
        int link_id = (vlinks_map.size() > 0) ? (int) vlinks_map.size() : 0;
        int link_src(stoi(it.first[0])), link_dest(stoi(it.first[1]));
        long link_bw(stol(it.second[0])), link_ltc(stol(it.second[1]));
        string link_desc("Link ");
        link_desc += myvnc->get_vnf_by_id(link_src)->description + '-' + myvnc->get_vnf_by_id(link_dest)->description;

        vlinks_map.insert({link_id, {link_id,link_desc,link_bw,link_ltc,myvnc->get_vnf_by_id(link_src),myvnc->get_vnf_by_id(link_dest)}});
    }

    cout <<  "\t====> [Request" << i << "] " << vlinks_map.size() << " VNC link(s) entries successfully loaded" << endl;

}
