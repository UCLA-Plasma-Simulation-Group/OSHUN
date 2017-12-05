/*!\brief  Parser routines - Definitions
* \author  PICKSC
 * \date   2017
 * \file   parser.cpp
 *
 * In here are the initialization routines for the State object.
 *
 */

//  Standard libraries
#include <iostream>
#include <vector>
#include <valarray>
#include <complex>
#include <algorithm>
#include <cstdlib>
#include <cfloat>

#include <math.h>
#include <map>

//  My libraries
#include "lib-array.h"
#include "lib-algorithms.h"
#include "external/exprtk.hpp"

//  Declerations
#include "parser.h"

//----------------------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------------------
void Parser::checkparse(parser_t& parser, std::string& expression_str, expression_t& expression){
//----------------------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------------------
    if (!parser.compile(expression_str,expression))
    {
        printf("Error: %s\tExpression: %s\n",
               parser.error().c_str(),
               expression_str.c_str());

        for (std::size_t i = 0; i < parser.error_count(); ++i)
        {
            error_t error = parser.get_error(i);
            printf("Error: %02d Position: %02d Type: [%s] Msg: %s Expr: %s\n",
                   static_cast<int>(i),
                   static_cast<int>(error.token.position),
                   exprtk::parser_error::to_str(error.mode).c_str(),
                   error.diagnostic.c_str(),
                   expression_str.c_str());
        }

        exit(1);
    }

}
//----------------------------------------------------------------------------------------------------------------------------
/**
 * @brief      parses function over x (or a grid) and returns a grid.
 *
 * @param      str_profile  The string profile
 * @param      profile      The profile
 */
//----------------------------------------------------------------------------------------------------------------------------
void Parser::parseprofile( const valarray<double>& grid, std::string& str_profile, valarray<double>& profile){

    /// Find curly brackets
    std::size_t posL = str_profile.find("{");
    std::size_t posR = str_profile.find("}");

    /// Uniform profile, convert string to double and set.
    if (str_profile[0] == 'c')
    {
        // std::cout<< "\nexpression_str is  " << str_profile.substr(posL+1,posR-(posL+1)) << "\n";
        
        // profile = std::stod(str_profile.substr(posL+1,posR-(posL+1)));
        profile = std::strtod((str_profile.substr(posL+1,posR-(posL+1))).c_str(),NULL);
        

        // std::cout << "\nprofile[1] = " << writeable[0] << " , " << writeable[1] << "\n\n";
        // std::cout << "\nprofile[1]= " << profile[1] << "\n\n";
    }

    // exit(1);

    /// Profile defined by function
    if (str_profile[0] == 'f')
    {
        /// Grab everything inside curly brackets
        std::string expression_str = str_profile.substr(posL+1,posR-(posL+1));

        symbol_table_t symbol_table;
        symbol_table.add_constants();
        double xstr;
        symbol_table.add_variable("x",xstr);

        expression_t expression;

        expression.register_symbol_table(symbol_table);

        parser_t parser;

        checkparse(parser, expression_str, expression);

        for (size_t i(0); i < profile.size(); ++i) {
            xstr = grid[i];
            profile[i] = expression.value();
        }
    }

    if (str_profile[0] == 'p'){
        std::string pcw_pairs = str_profile.substr(posL+1,posR-(posL+1));
        vector<size_t> positions_comma, positions_semicolon; // holds all the positions that sub occurs within str
        vector<double> loc,val;
        size_t pos = pcw_pairs.find(',', 0);
        while(pos != string::npos)
        {
            positions_comma.push_back(pos);
            pos = pcw_pairs.find(',',pos+1);
        }

        pos = pcw_pairs.find(';', 0);
        while(pos != string::npos)
        {
            positions_semicolon.push_back(pos);
            pos = pcw_pairs.find(';',pos+1);
        }

        if (positions_comma.size() != positions_semicolon.size()){
            std::cout << "Error in piecewise input, the values must be supplied in pairs with ',' in between the two values in each pair, and ending each pair with a ';'. E.g. {0.1,1.0;0.2,1.5;} \n";

            exit(1);
        }


        // loc.push_back(std::stod(pcw_pairs.substr(0,positions_comma[0]-1)));
        // std::string
        loc.push_back(std::strtod((pcw_pairs.substr(0,positions_comma[0]-1)).c_str(),NULL));
        // val.push_back(std::stod(pcw_pairs.substr(positions_comma[0]+1,positions_semicolon[0]-positions_comma[0]-1)));
        val.push_back(std::strtod((pcw_pairs.substr(positions_comma[0]+1,positions_semicolon[0]-positions_comma[0]-1)).c_str(),NULL));

        // std::cout << "loc0 = " << loc[0] << endl;
        // std::cout << "val0 = " << val[0] << endl;

        for (size_t i(1); i < positions_comma.size(); ++i) {

            // loc.push_back(std::stod(pcw_pairs.substr(positions_semicolon[i-1]+1,positions_comma[i]-positions_semicolon[i-1]-1)));
            loc.push_back(std::strtod((pcw_pairs.substr(positions_semicolon[i-1]+1,positions_comma[i]-positions_semicolon[i-1]-1)).c_str(),NULL));

            // val.push_back(std::stod(pcw_pairs.substr(positions_comma[i]+1,positions_semicolon[i]-positions_comma[i]-1)));
            val.push_back(std::strtod((pcw_pairs.substr(positions_comma[i]+1,positions_semicolon[i]-positions_comma[i]-1)).c_str(),NULL));

            // std::cout << "loci = " << loc[i] << endl;
            // std::cout << "vali = " << val[i] << endl;

        }


        double tmp_x, xx;
        size_t t(0);
        for (size_t i(0); i < grid.size(); ++i) {
            t = 1;
            xx = grid[i];
            if (grid[i] < loc[0])               xx += loc[loc.size()-1] - loc[0];
            if (grid[i] > loc[loc.size()-1])  xx += loc[0]-loc[loc.size()-1] ;

            while ((xx > loc[t]) && (t < loc.size()-1)) ++t;
            tmp_x = val[t-1] + (val[t]-val[t-1])/(loc[t]-loc[t-1])
                               * (xx - loc[t-1]);
            // std::cout << i << "     " << x[i] << "      " << tmp_x << std::endl;
            //             for (size_t j(0); j < Temp.dim2(); ++j)  Temp(i,j) *= tmp_x;
            profile[i] = tmp_x;
        }
    }
}

//----------------------------------------------------------------------------------------------------------------------------
/**
 * @brief      parses function over x (or a grid) and returns a grid.
 *
 * @param      str_profile  The string profile
 * @param      profile      The profile
 */
//----------------------------------------------------------------------------------------------------------------------------
void Parser::parseprofile( const valarray<double>& grid1, const valarray<double>& grid2, std::string& str_profile, Array2D<double>& profile){

    /// Find curly brackets
    std::size_t posL = str_profile.find("{");
    std::size_t posR = str_profile.find("}");

    /// Uniform profile, convert string to double and set.
    if (str_profile[0] == 'c')
    {
        // std::cout<< "\nexpression_str is  " << str_profile.substr(posL+1,posR-(posL+1)) << "\n";
        
        // profile = std::stod(str_profile.substr(posL+1,posR-(posL+1)));
        profile = std::strtod((str_profile.substr(posL+1,posR-(posL+1))).c_str(),NULL);
        

        // std::cout << "\nprofile[1] = " << writeable[0] << " , " << writeable[1] << "\n\n";
        // std::cout << "\nprofile[1]= " << profile[1] << "\n\n";
    }

    /// Profile defined by function
    if (str_profile[0] == 'f')
    {
        /// Grab everything inside curly brackets
        std::string expression_str = str_profile.substr(posL+1,posR-(posL+1));

        symbol_table_t symbol_table;
        symbol_table.add_constants();
        double xstr, ystr;
        symbol_table.add_variable("x",xstr);
        symbol_table.add_variable("y",ystr);

        expression_t expression;

        expression.register_symbol_table(symbol_table);

        parser_t parser;

        checkparse(parser, expression_str, expression);

        for (size_t i1(0); i1 < grid1.size(); ++i1) 
        {
            xstr = grid1[i1];
            
            for (size_t i2(0); i2 < grid2.size(); ++i2) 
            {
                ystr = grid2[i2];
                
                profile(i1,i2) = expression.value();

            }
        }
    }

    // if (str_profile[0] == 'p')
    // {
    //     std::string pcw_pairs = str_profile.substr(posL+1,posR-(posL+1));
    //     vector<size_t> positions_comma, positions_semicolon; // holds all the positions that sub occurs within str
    //     vector<double> loc,val;
    //     size_t pos = pcw_pairs.find(',', 0);
    //     while(pos != string::npos)
    //     {
    //         positions_comma.push_back(pos);
    //         pos = pcw_pairs.find(',',pos+1);
    //     }

    //     pos = pcw_pairs.find(';', 0);
    //     while(pos != string::npos)
    //     {
    //         positions_semicolon.push_back(pos);
    //         pos = pcw_pairs.find(';',pos+1);
    //     }

    //     if (positions_comma.size() != positions_semicolon.size()){
    //         std::cout << "Error in piecewise input, the values must be supplied in pairs with ',' in between the two values in each pair, and ending each pair with a ';'. E.g. {0.1,1.0;0.2,1.5;} \n";

    //         exit(1);
    //     }


    //     // loc.push_back(std::stod(pcw_pairs.substr(0,positions_comma[0]-1)));
    //     // std::string
    //     loc.push_back(std::strtod((pcw_pairs.substr(0,positions_comma[0]-1)).c_str(),NULL));
    //     // val.push_back(std::stod(pcw_pairs.substr(positions_comma[0]+1,positions_semicolon[0]-positions_comma[0]-1)));
    //     val.push_back(std::strtod((pcw_pairs.substr(positions_comma[0]+1,positions_semicolon[0]-positions_comma[0]-1)).c_str(),NULL));

    //     // std::cout << "loc0 = " << loc[0] << endl;
    //     // std::cout << "val0 = " << val[0] << endl;

    //     for (size_t i(1); i < positions_comma.size(); ++i) {

    //         // loc.push_back(std::stod(pcw_pairs.substr(positions_semicolon[i-1]+1,positions_comma[i]-positions_semicolon[i-1]-1)));
    //         loc.push_back(std::strtod((pcw_pairs.substr(positions_semicolon[i-1]+1,positions_comma[i]-positions_semicolon[i-1]-1)).c_str(),NULL));

    //         // val.push_back(std::stod(pcw_pairs.substr(positions_comma[i]+1,positions_semicolon[i]-positions_comma[i]-1)));
    //         val.push_back(std::strtod((pcw_pairs.substr(positions_comma[i]+1,positions_semicolon[i]-positions_comma[i]-1)).c_str(),NULL));

    //         // std::cout << "loci = " << loc[i] << endl;
    //         // std::cout << "vali = " << val[i] << endl;

    //     }


    //     double tmp_x, xx;
    //     size_t t(0);
    //     for (size_t i(0); i < grid.size(); ++i) {
    //         t = 1;
    //         xx = grid[i];
    //         if (grid[i] < loc[0])               xx += loc[loc.size()-1] - loc[0];
    //         if (grid[i] > loc[loc.size()-1])  xx += loc[0]-loc[loc.size()-1] ;

    //         while ((xx > loc[t]) && (t < loc.size()-1)) ++t;
    //         tmp_x = val[t-1] + (val[t]-val[t-1])/(loc[t]-loc[t-1])
    //                            * (xx - loc[t-1]);
    //         // std::cout << i << "     " << x[i] << "      " << tmp_x << std::endl;
    //         //             for (size_t j(0); j < Temp.dim2(); ++j)  Temp(i,j) *= tmp_x;
    //         profile[i] = tmp_x;
    //     }
    // }
}
//----------------------------------------------------------------------------------------------------------------------------
/**
 * @brief      parses function over t (or single variable) and returns a value
 */
//----------------------------------------------------------------------------------------------------------------------------
void Parser::parseprofile(const double& input, std::string& str_profile, double& output){

    /// Find curly brackets
    std::size_t posL = str_profile.find("{");
    std::size_t posR = str_profile.find("}");

    /// Uniform profile, convert string to double and set.
    if (str_profile[0] == 'c')
    {
        // std::cout<< "expression_str is  " << str_profile.substr(posL+1,posR-(posL+1)) << "\n";
        // output = std::stod(str_profile.substr(posL+1,posR-(posL+1)));
        output = std::strtod((str_profile.substr(posL+1,posR-(posL+1))).c_str(),NULL);
        
    }

    /// Profile defined by function
    if (str_profile[0] == 'f')
    {
        /// Grab everything inside curly brackets
        std::string expression_str = str_profile.substr(posL+1,posR-(posL+1));

        symbol_table_t symbol_table;
        symbol_table.add_constants();
        double xstr;
        symbol_table.add_variable("t",xstr);

        expression_t expression;

        expression.register_symbol_table(symbol_table);

        parser_t parser;

        checkparse(parser, expression_str, expression);

        // for (size_t i(0); i < profile.size(); ++i) {
        xstr = input;
        // std::cout<< "i= " << i << ",xstr = " << xstr << "\n";

        // result =
        // Temperature_map(s,i) = expression_Temperature.value();
        output = expression.value();
        // std::cout  << " \n 10: " << input << " , " << output << "\n";
        // exit(1);
        // printf("Result: %10.5f\n",expression.value());

    }

    //     if (str_profile[0] == 'p'){
    //         std::string pcw_pairs = str_profile.substr(posL+1,posR-(posL+1));
    //         vector<size_t> positions_comma, positions_semicolon; // holds all the positions that sub occurs within str
    //         vector<double> loc,val;
    //         size_t pos = pcw_pairs.find(',', 0);
    //         while(pos != string::npos)
    //         {
    //             positions_comma.push_back(pos);
    //             pos = pcw_pairs.find(',',pos+1);
    //         }

    //         pos = pcw_pairs.find(';', 0);
    //         while(pos != string::npos)
    //         {
    //             positions_semicolon.push_back(pos);
    //             pos = pcw_pairs.find(';',pos+1);
    //         }

    //         if (positions_comma.size() != positions_semicolon.size()){
    //             std::cout << "Error in piecewise input, the values must be supplied in pairs with ',' in between the two values in each pair, and ending each pair with a ';'. E.g. {0.1,1.0;0.2,1.5;} \n";

    //             exit(1);
    //         }


    //         loc.push_back(std::stod(pcw_pairs.substr(0,positions_comma[0]-1)));
    //         val.push_back(std::stod(pcw_pairs.substr(positions_comma[0]+1,positions_semicolon[0]-positions_comma[0]-1)));

    //         // std::cout << "loc0 = " << loc[0] << endl;
    //         // std::cout << "val0 = " << val[0] << endl;

    //         for (size_t i(1); i < positions_comma.size(); ++i) {

    //                 loc.push_back(std::stod(pcw_pairs.substr(positions_semicolon[i-1]+1,positions_comma[i]-positions_semicolon[i-1]-1)));
    //                 val.push_back(std::stod(pcw_pairs.substr(positions_comma[i]+1,positions_semicolon[i]-positions_comma[i]-1)));

    //                 // std::cout << "loci = " << loc[i] << endl;
    //                 // std::cout << "vali = " << val[i] << endl;

    //         }


    //         double tmp_x, xx;
    //         size_t t(0);
    //         for (size_t i(0); i < grid.size(); ++i) {
    //             t = 1;
    //             xx = grid[i];
    //             if (grid[i] < loc[0])               xx += loc[loc.size()-1] - loc[0];
    //             if (grid[i] > loc[loc.size()-1])  xx += loc[0]-loc[loc.size()-1] ;

    //             while ((xx > loc[t]) && (t < loc.size()-1)) ++t;
    //             tmp_x = val[t-1] + (val[t]-val[t-1])/(loc[t]-loc[t-1])
    //                                 * (xx - loc[t-1]);
    //             // std::cout << i << "     " << x[i] << "      " << tmp_x << std::endl;
    // //             for (size_t j(0); j < Temp.dim2(); ++j)  Temp(i,j) *= tmp_x;
    //             profile[i] = tmp_x;
    //         }
    //     }
}
//----------------------------------------------------------------------------------------------------------------------------
/**
 * @brief      parses a grid and a variable function (x and t) and returns a grid that's modulated in time
 *
 */
//----------------------------------------------------------------------------------------------------------------------------
void Parser::parseprofile(const valarray<double>& grid, const double& input, std::string& str_profile, valarray<double>& output){

    /// Find curly brackets
    std::size_t posL = str_profile.find("{");
    std::size_t posR = str_profile.find("}");

    /// Uniform profile, convert string to double and set.
    if (str_profile[0] == 'c')
    {
        // std::cout<< "expression_str is  " << str_profile.substr(posL+1,posR-(posL+1)) << "\n";
        // output = std::stod(str_profile.substr(posL+1,posR-(posL+1)));
        output = std::strtod((str_profile.substr(posL+1,posR-(posL+1))).c_str(),NULL);
    }

    /// Profile defined by function
    if (str_profile[0] == 'f')
    {
        /// Grab everything inside curly brackets
        std::string expression_str = str_profile.substr(posL+1,posR-(posL+1));

        symbol_table_t symbol_table;
        symbol_table.add_constants();
        double xstr, tstr;
        symbol_table.add_variable("x",xstr);
        symbol_table.add_variable("t",tstr);

        expression_t expression;

        expression.register_symbol_table(symbol_table);

        parser_t parser;

        checkparse(parser, expression_str, expression);
        tstr = input;

        for (size_t i(0); i < output.size(); ++i) {
            xstr = grid[i];

            output[i] = expression.value();
        }



    }

}

//----------------------------------------------------------------------------------------------------------------------------
void Parser::parseprofile(const valarray<double>& grid1, const valarray<double>& grid2, const double& input, std::string& str_profile, Array2D<double>& output){

    /// Find curly brackets
    std::size_t posL = str_profile.find("{");
    std::size_t posR = str_profile.find("}");

    /// Uniform profile, convert string to double and set.
    if (str_profile[0] == 'c')
    {
        // std::cout<< "expression_str is  " << str_profile.substr(posL+1,posR-(posL+1)) << "\n";
        // output = std::stod(str_profile.substr(posL+1,posR-(posL+1)));
        output = std::strtod((str_profile.substr(posL+1,posR-(posL+1))).c_str(),NULL);
    }

    /// Profile defined by function
    if (str_profile[0] == 'f')
    {
        /// Grab everything inside curly brackets
        std::string expression_str = str_profile.substr(posL+1,posR-(posL+1));

        symbol_table_t symbol_table;
        symbol_table.add_constants();
        double xstr, ystr, tstr;
        symbol_table.add_variable("x",xstr);
        symbol_table.add_variable("y",ystr);
        symbol_table.add_variable("t",tstr);

        expression_t expression;

        expression.register_symbol_table(symbol_table);

        parser_t parser;

        checkparse(parser, expression_str, expression);
        tstr = input;

        for (size_t ix(0); ix < output.dim1(); ++ix) 
        {
            for (size_t iy(0); iy < output.dim2(); ++iy) 
            {

                xstr = grid1[ix];
                ystr = grid2[iy];

                output(ix,iy) = expression.value();
            }
        }



    }

}

