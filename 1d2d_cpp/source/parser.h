
/*! \brief Function Parser - Declarations
 * \author PICKSC
 * \date   2017
 * \file   input.h
 * 
 * Contains:
 * 1) input reader
 * 2) default values for input variables.
 */

#ifndef DECL_PARSER_H
#define DECL_PARSER_H



//**************************************************************
//--------------------------------------------------------------
namespace Parser{
    typedef exprtk::symbol_table<double> symbol_table_t;
    typedef exprtk::expression<double> expression_t;
    typedef exprtk::parser<double> parser_t;
    typedef exprtk::parser_error::type error_t;

    void checkparse(parser_t& parser, std::string& expression_str, expression_t& expression);
    
    void parseprofile(const std::valarray<double>& grid, std::string& str_profile, std::valarray<double>& profile);
    void parseprofile(const std::valarray<double>& gridx, const std::valarray<double>& gridy, std::string& str_profile, Array2D<double>& profile);

    void parseprofile(const double& input, std::string& str_profile, double& ouput);
	void parseprofile(const std::valarray<double>& grid, const double& input, std::string& str_profile, std::valarray<double>& ouput);
	void parseprofile(const std::valarray<double>& grid1, const std::valarray<double>& grid2, const double& input, std::string& str_profile, Array2D<double>& ouput);
}

#endif
