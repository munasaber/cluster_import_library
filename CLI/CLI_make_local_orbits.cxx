#include <CLI/CLI.hpp>
#include "CLI/App.hpp"


int main(int argc, char* argv[]) { 
	//Provides app description when the -h tag is applied
	CLI::App app{"This allows one to find cluster coordinates and indicies of a particular structure based on cutoff radius"};
	
	//Tag for the user input for the structures path
	casmutils::fs::path structurepath;
	app.add_option("-s, --structure", structurepath, "Please input the file path of the base structure")-> required();

	//Tag for the user input of .json file path
	casmutils::fs::path jsonpath;
	app.add_option("-j, --json_path", jsonpath, "Please input the path to the input json file");

	//Tag for the path to the output file (This is not required so otherwise the file will be output to the screen)
	casmutils::fs::path outputpath
	app.add_option("-o, --output", outputpath, "Please input the path to the output file. If this option is not used, the cluster data will be printed to the screen"); 

}
