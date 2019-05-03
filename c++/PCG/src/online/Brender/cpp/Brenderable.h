/*
* @author: Gustavo Lopez 10-21-17
*
* @version: 1.0
*/

#pragma once
#include <string>
#include <fstream>
#include <vector>
#include <memory>
#include "BrenderManager.h"


class Brenderable
{
public:
	enum BranderableType {
		Truncate, Append, ResetAppend
	};

	Brenderable() {};
	virtual ~Brenderable() {};
	virtual int getBrenderCount() const { return 1; };
	virtual std::vector<std::string> getBrenderNames() const { return std::vector<std::string>(1, ""); };
	virtual std::vector<std::string> getBrenderExtensions() const { return std::vector<std::string>(1, ""); };
	virtual std::vector<int> getBrenderTypes() const { return std::vector<int>(1, Truncate); };
	virtual void exportBrender(std::vector< std::shared_ptr< std::ofstream > > outfiles, int frame, double time) const = 0;
private:

};