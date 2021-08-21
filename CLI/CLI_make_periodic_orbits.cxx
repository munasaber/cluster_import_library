#include <CLI/Formatter.hpp>
#include <CLI/Config.hpp>
#include <CLI/CLI.hpp>
#include <CLI/App.hpp>
#include <casm/casm_io/json/jsonParser.hh>
#include <casm/casm_io/json/jsonFile.hh>
#include <casmutils/definitions.hpp>
#include <casmutils/xtal/structure.hpp>
#include <casmutils/xtal/symmetry.hpp>
#include <casmutils/sym/cartesian.hpp>
#include <casmutils/xtal/site.hpp>
#include <casmutils/xtal/structure_tools.hpp>
#include <casmutils/clusterography/cluster_finder.hpp>
#include <vector>
#include <nlohmann/json.hpp>
#include <../tests/autotools.hh>
#include "../../CASMcode/include/casm/external/Eigen/Core"
#include "../../CASMcode/include/casm/external/Eigen/Dense"
using json= nlohmann::json;


//Gets the labels of the orbit sizes that you want to find the clusters for from an input json file
std::vector<int> get_orbit_labels(casmutils::fs::path& jsonpath)
{
	std::ifstream ifs(jsonpath);
	json j = json::parse(ifs);
	return j.at("orbit_branch_specs");
}


//get maximum length of each radius from an input json file
std::vector<double> get_max_lengths(casmutils::fs::path& jsonpath)
{
	std::ifstream ifs(jsonpath);
	json j = json::parse(ifs);
	return j.at("max_lengths");
}

//necessary for inputting the filtered_sites
bool isDivisibleby4(int value){
   if(value % 4 == 0)
   {
      return true;
   }
   return false;
}


//get the sites that the user wants to filter out as stated in the site_filter label in the input json file
std::vector<casmutils::xtal::Site> get_filtered_sites(casmutils::fs::path& jsonpath)
{
	std::ifstream ifs(jsonpath);
	json j = json::parse(ifs);
	//std::vector<double> vector_of_sites=j.at("site_filter");
	std::vector<std::string> vector_of_sites=j.at("site_filter");
	//std::vector<Eigen::Vector3d> total_vector;
	Eigen::Vector3d my_eigen;
	std::string occupant_type;
	std::vector<casmutils::xtal::Site> Site_list;
	for (int i=0; i< vector_of_sites.size(); i=i+4)
	{	
		Eigen::Vector3d my_eigen(std::stod(vector_of_sites[i]), std::stod(vector_of_sites[i+1]), std::stod(vector_of_sites[i+2]));
		Site_list.emplace_back(casmutils::xtal::Site(my_eigen, vector_of_sites[i+3]));
	}
	if( !isDivisibleby4(vector_of_sites.size()) ) 
	{
      		throw "Division by zero condition!";
        }
	return Site_list;
}

//make function to use as input to cluster function
std::function<bool(const casmutils::xtal::Site)> make_site_filter_func(std::vector<casmutils::xtal::Site> sites_to_be_excluded)
{
	std::function<bool(const casmutils::xtal::Site)> my_func=[&](casmutils::xtal::Site site)
	{
		casmutils::xtal::SiteEquals_f site_comparator(1e-6);
		for (const auto&  x: sites_to_be_excluded)
		{
			if (site_comparator(x, site))
			{
				return false;
			}
		}
	return true;
	};
return my_func;
}


//Convert json to data and manipulate using the cluster_finder functions
std::vector<casmutils::cluster::Orbit> make_orbits(casmutils::fs::path& jsonpath, casmutils::fs::path& structurepath)
{
	
	//Get maxlength from json
	const std::vector<double> maxlengths=get_max_lengths(jsonpath);

	//get structure pointer from structure path
	
	casmutils::xtal::Structure structure= casmutils::xtal::Structure::from_poscar(structurepath);
	
         	
	//get left out clusters from user input
	const std::vector<casmutils::xtal::Site> filtered_sites(get_filtered_sites(jsonpath));
	std::vector<casmutils::cluster::Orbit> total_orbits= casmutils::cluster::make_periodic_orbits(maxlengths, structure, make_site_filter_func(filtered_sites));
	return total_orbits;
	
};






int main(int argc, char* argv[]) { 
	
	//Provides app description when the -h tag is applied
	CLI::App app{"This allows one to find cluster coordinates and indicies of a particular structure based on cutoff radius"};
	
	//Tag for the user input for the structures path
        casmutils::fs::path structurepath;
	CLI::Option* structure_option=app.add_option("-s, --structure", structurepath, "Please input the file path of the base structure");
        
	//Tag for maximum cutoff radius for sites to be included into the cluster
	std::vector<double>  maxlengths;
	CLI::Option* length_option=app.add_option("-m, --max_length", maxlengths, "Select the maximum cutoff radius for sites to be included in a given cluster. Higher order clusters must either have a smaller or same size radius as lower order clusters.");
	
	//Tag for the user input of .json file path
	casmutils::fs::path jsonpath;
	CLI::Option* json_option=app.add_option("-j, --json_path", jsonpath, "Please input the path to the input json file");
		
	//Tag for the path to the output file (This is not required so otherwise the file will be output to the screen)
	casmutils::fs::path outputpath;
	CLI::Option* output_option=app.add_option("-o, --output", outputpath, "Please input the path to the output file. If this option is not used, the cluster data will be printed to the screen"); 
        CLI::Option* desc_option=app.add_flag("--desc", "Prints to the screen an explanation for te json format"); 

	CLI11_PARSE(app, argc, argv);
	
	if (* desc_option)
	{
		std::cout<<"Input json requires an 'orbit_branch_specs' input that gives the size of each cluster type (i.e. branch)"<<std::endl;
		std::cout<<"Clusters listed in 'orbit_branch_specs' must be consecutive starting from 0 as each branch is dependent on the last"<<std::endl;
		std::cout<<"A maximum length for each cluster type must be listed with the flag 'max_length'"<<std::endl;
		std::cout<<"Optionally, selected sites can be ignored in the determination of orbits. To allow for this, use the 'site_filter' tag"<<std::endl;
		std::cout<<"A sample json file for a cluster up to quadruplet is given below: "<<std::endl;
		std::cout<<"{"<<std::endl;
		std::cout<<"     \"orbit_branch_specs\" : [0, 1, 2, 3, 4],"<<std::endl;
		std::cout<<"     \"max_lengths\" : [5, 5, 5, 5, 5],"<<std::endl;	
		std::cout<<"     \"site_filter\" : [\"0.00000000\", \"0.00000000\", \"0.00000000\", \"A\"]"<<std::endl;
		std::cout<<"}"<<std::endl;
		return 0;	
	}
	std::vector<int> orbit_labels=get_orbit_labels(jsonpath);
	std::vector<double> max_lengths=get_max_lengths(jsonpath);
	std::vector<casmutils::cluster::Orbit> all_filtered_orbits= make_orbits(jsonpath, structurepath);
	//If the user inputs a json, this loop prints out the data to the users outpathfile path of choice
	if (* output_option)
	{
		std::ofstream outpathfile(outputpath);
		outpathfile<<"Here are the coordinates for each set of clusters given the specified lengths, coordinate dimensions, and filtered site information"<<std::endl;
		int b=0;
		int i=0;
		for (const auto& orbit: all_filtered_orbits)
		{
			
			outpathfile<<"    "<<i<<" of "<<all_filtered_orbits.size()<<" Orbits     Orbit size: "<<orbit.size()<<std::endl;
			std::vector<casmutils::cluster::Cluster> clusters=orbit.get_clusters();
			int j=0;
			for (const auto& cluster : clusters)
			{
				outpathfile<<"        Sites in "<<j<<" of "<<cluster.size()<<" Equivalent Clusters in Orbit "<<i<<std::endl;
				std::vector<casmutils::xtal::Site> sites=cluster.get_sites();
				for (const auto& site: sites)
				{
				  outpathfile<<"             "<<site.frac(casmutils::xtal::Structure::from_poscar(structurepath).lattice()).transpose()<<" "<<site.label()<<std::endl;
				}
				j++;
				//Gets each Site coordinate from each Cluster in Each orbit
			}
			outpathfile<<std::endl;
			i++;
		}
	}
	//If output option is not selected then print Sites of CLusterrs of Orbits to the screen
	else 
	{

		std::ofstream outpathfile(outputpath);
		std::cout<<"Here are the coordinates for each set of clusters given the specified lengths, coordinate dimensions, and filtered site information"<<std::endl;
		int i=0;
		for (const auto& orbit: all_filtered_orbits)
		{
			std::cout<<i<<" of "<<all_filtered_orbits.size()<<" Orbits     Orbit size: "<<orbit.size()<<std::endl;
			std::vector<casmutils::cluster::Cluster> clusters=orbit.get_clusters();
			int j=0;
			for (const auto& cluster : clusters)
			{
				std::cout<<"        Sites in "<<j<<" of "<<cluster.size()<<" Equivalent Clusters in Orbit "<<i<<std::endl;
				std::vector<casmutils::xtal::Site> sites=cluster.get_sites();
				for (const auto& site: sites)
				{
					std::cout<<"             "<<site.frac(casmutils::xtal::Structure::from_poscar(structurepath).lattice()).transpose()<<" "<<site.label()<<std::endl;
				}
				j++;
				//Gets each Site coordinate from each Cluster in Each orbit
			}
			std::cout<<std::endl;
			i++;
		}


	}
	
	
	return 0; 
}
