#include "Parser.h"

Parser::Parser(std::string filename):m_filename(filename) 
{
    std::cout<<"\t\tREADINE FILE "<<m_filename<<"\n";

    auto cst_ZERO = [](double u, double v, double t){
        return 0.0;
    };

    for (size_t idx=0; idx<5; idx++)
    {
        param.Functions.push_back(cst_ZERO);
    }

}


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
            currentSection = trimBlank(cleanLine.substr(1, cleanLine.size() - 2));
            std::transform(currentSection.begin(), currentSection.end(), currentSection.begin(), ::toupper);
        }

        size_t eq_pos = cleanLine.find('=');
        if (eq_pos == std::string::npos) continue;

        std::string key = trimBlank(cleanLine.substr(0, eq_pos));
        std::transform(key.begin(), key.end(), key.begin(), ::toupper);
        std::string value = trimBlank(cleanLine.substr(eq_pos+1));
        std::transform(value.begin(), value.end(), value.begin(), ::toupper);
        if (currentSection=="SOLVER") m_rawSolverParameters[key]=value;
        else if (currentSection=="GRID") m_rawGridParameters[key]=value;
        else if(currentSection=="BOUNDARYCONDITIONS") m_rawBoundaryConditions[key]=value;
        else if (currentSection=="FUNCTIONS") m_rawFunctions[key]=value;
    }
}

void Parser::convertToParameters()
{
    

    try{
        param.x_start = std::stod(m_rawGridParameters.at("X_START"));
        param.x_end = std::stod(m_rawGridParameters.at("X_END"));
        param.growFactor = std::stod(m_rawGridParameters.at("GROW_FACTOR"));
        param.A = std::stod(m_rawGridParameters.at("A"));
        param.integrationOrderGrid = std::stoi(m_rawGridParameters.at("INTEGRATION_ORDER"));
        param.numbOfElements = std::stoi(m_rawGridParameters.at("NUMB_OF_ELEMENTS"));

    }
    catch (const std::exception &e)
    {
        std::cout<<"[ ERROR IN GRID SECTION ] "<<e.what()<<"\n";
    }

    try{
        param.work_type = m_rawSolverParameters.at("WORK_TYPE");
        param.beta = std::stod(m_rawSolverParameters.at("BETA"));
        param.gamma = std::stod(m_rawSolverParameters.at("GAMMA"));
        param.epsilon = std::stod(m_rawSolverParameters.at("EPSILON"));
        param.omega = std::stod(m_rawSolverParameters.at("OMEGA"));
        param.t_max = std::stod(m_rawSolverParameters.at("T_MAX"));
    }
    catch (const std::exception &e)
    {
        std::cout<<"[ ERROR IN SOLVER SECTION ] "<<e.what()<<"\n";
    }

    try{
        for (const auto& [key, value] : m_rawBoundaryConditions)
        {
            if(key.substr(0,2)=="BC")
            {
                std::stringstream ss(value);
                std::string segment;
                std::vector<std::string> parts;

                while(std::getline(ss, segment, '/')) 
                {
                    parts.push_back(trimBlank(segment));
                }

                if (parts.size()!=3) continue;
                int boundaryConditionType;
                
                if(parts[0]=="DIRICHLET") boundaryConditionType=1;
                else if (parts[0]=="NEUMAN") boundaryConditionType=2;

                param.boundaryConditions.push_back({std::stoi(parts[1]), boundaryConditionType, std::stod(parts[2])});
            }
        }
    }
    catch (const std::exception &e)
    {
        std::cout<<"[ ERROR IN BOUNDARY CONDITIONS SECTION ] "<<e.what()<<"\n";
    }

    try{

        struct EquationContext{
            double x=0.0;
            double u=0.0;
            double t=0.0;
            exprtk::symbol_table<double> symbol_table;
            exprtk::expression<double> expression;

        };
        for (const auto& [key,value] : m_rawFunctions)
        {
            auto ctx = std::make_shared<EquationContext>();

            ctx->symbol_table.add_variable("X", ctx->x);
            ctx->symbol_table.add_variable("U", ctx->u);
            ctx->symbol_table.add_variable("T", ctx->t);
            ctx->symbol_table.add_constants();

            ctx->expression.register_symbol_table(ctx->symbol_table);
            exprtk::parser<double> expr_parser;

            if (!expr_parser.compile(value,ctx->expression)) throw std::runtime_error("Compilation error "+value);

            auto func = [ctx](double x_val, double u_val, double t_val) -> double
            {
                ctx->x=x_val;
                ctx->u=u_val;
                ctx->t=t_val;

                return ctx->expression.value();
            };

            if (key=="B")
                param.Functions[0]=func;
            else if (key=="C")
                param.Functions[1]=func;
            else if (key=="D")
                param.Functions[2]=func;
            else if (key=="E")
                param.Functions[3]=func;
            else if (key=="F")
                param.Functions[4]=func;
            else
                throw std::runtime_error("UNKNOWN FUNCTION SYMBOL: "+key);


        }

    }
    catch (const std::exception &e)
    {
        std::cout<<"[ ERROR IN FUNCTIONS SECTION ]"<<e.what()<<"\n";
    }
}

ConfigParameters Parser::getParameters() {return param;}