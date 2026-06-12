#include "Parser.h"

Parser::Parser(std::string filename):m_filename(filename) {std::cout<<"\t\tREADINE FILE "<<m_filename<<"\n";}


std::string Parser::trimBlank(const std::string &str)
{
    size_t first = str.find_first_not_of(" \t\r\n");
    size_t last = str.find_last_not_of(" \t\r\n");
    return (first == std::string::npos) ? "" : str.substr(first, (last - first + 1));
}

void Parser::readRawData()
{
    std::ifstream file(m_filename);
    std::string line;
    std::string currentSection;

    while(std::getline(file, line))
    {
        std::string cleanLine=trimBlank(line);
        if (cleanLine[0]=='#' || cleanLine.empty()) continue;

        if (cleanLine.front() == '[' && cleanLine.back() == ']') 
        {
            currentSection = cleanLine.substr(1, cleanLine.size() - 2);
            std::transform(currentSection.begin(), currentSection.end(), currentSection.begin(), ::toupper);
        }

        size_t eq_pos = line.find('=');
        if (eq_pos=std::string::npos) continue;

        std::string key = line.substr(0, eq_pos);
        std::transform(key.begin(), key.end(), key.begin(), ::toupper);
        std::string value = line.substr(eq_pos+1);
        std::transform(value.begin(), value.end(), value.begin(), ::toupper);
        if (currentSection=="SOLVER") m_rawSolverParameters[key]=value;
        else if (currentSection=="GRID") m_rawGridParameters[key]=value;
        else if(currentSection=="BOUNDARYCONDITIONS") m_rawBoundaryConditions[key]=value;
    }
}

void Parser::convertToParameters()
{
    ConfigParameters param;

    try{
        param.x_start = std::stod(m_rawGridParameters.at("X_START"));
        param.x_end = std::stod(m_rawGridParameters.at("X_END"));
        param.growFactor = std::stod(m_rawGridParameters.at("GROW_FACTOR"));
        param.integrationOrderGrid = std::stod(m_rawGridParameters.at("INTEGRATION_ORDER"));
        param.numbOfElements = std::stod(m_rawGridParameters.at("NUMB_OF_ELEMENTS"));

    }
    catch (const std::exception &e)
    {
        std::cout<<"[ ERROR ]"<<e.what()<<"\n";
    }
}