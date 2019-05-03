/*
* @author: Gustavo Lopez 10-21-17
*
* @version: 1.0
*/

#include "BrenderManager.h"
#include "Brenderable.h"


using namespace std;

bool BrenderManager::instanceFlag_ = false;
BrenderManager* BrenderManager::manager_ = NULL;
BrenderManager* BrenderManager::getInstance()
{
	if (!instanceFlag_)
	{
		manager_ = new BrenderManager();
		instanceFlag_ = true;
		return manager_;
	}
	else
	{
		return manager_;
	}
}

int BrenderManager::getFrame() const
{
	return frame_;
}

void BrenderManager::exportBrender(double time)
{
	int objNum = 1;
	for (auto brenderable : brenderables_) {
		vector<string> names = brenderable->getBrenderNames();
		vector<string> extensions = brenderable->getBrenderExtensions();
		vector<int> types = brenderable->getBrenderTypes();
		vector< shared_ptr< ofstream > > outfiles;
		// Initialize files
		for (int i = 0; i < brenderable->getBrenderCount(); ++i) {
			auto outfile = make_shared<ofstream>();
			outfiles.push_back(outfile);

			char filename[512];

			if (extensions[i].compare("") == 0) extensions[i] = "obj";

			if (types[i] == Brenderable::Truncate) {
				// if object has not been given name
				if (names[i].compare("") == 0) {
					sprintf(filename, "%s/%06d_Object%d.%s", EXPORT_DIR_.c_str(), frame_, objNum, extensions[i].c_str());
				}
				// if object has been given specific name
				else {
					sprintf(filename, "%s/%06d_%s.%s", EXPORT_DIR_.c_str(), frame_, names[i].c_str(), extensions[i].c_str());
				}

				// open file
				outfile->open(filename, ofstream::out | ofstream::trunc);

				// frame string
				char framestring[50];
				sprintf(framestring, "# frame %06d \n", frame_);
				*outfile << framestring;
				// frame time
				char timeval[50];
				sprintf(timeval, "# time %f \n", time);
				*outfile << timeval;
				// obj name
				// if object has not been given name
				if (names[i].compare("") == 0) {
					*outfile << "# name Object" + to_string(objNum) + " \n";
				}
				// if object has been given specific name
				else {
					*outfile << "# name " + names[i] + " \n";
				}
			}
			else if (types[i] == Brenderable::Append || types[i] == Brenderable::ResetAppend) {
				// if object has not been given name
				if (names[i].compare("") == 0) {
					sprintf(filename, "%s/Object%d.%s", EXPORT_DIR_.c_str(), objNum, extensions[i].c_str());
				}
				// if object has been given specific name
				else {
					sprintf(filename, "%s/%s.%s", EXPORT_DIR_.c_str(), names[i].c_str(), extensions[i].c_str());
				}

				if (types[i] == Brenderable::Append) {
					outfile->open(filename, ofstream::out | ofstream::app);
				}
				else {
					outfile->open(filename, ofstream::out | ofstream::trunc);
				}
			}
		}
		// Write to files
		brenderable->exportBrender(outfiles, frame_, time);
		// Close files
		for (int i = 0; i < brenderable->getBrenderCount(); ++i) {
			outfiles[i]->close();
		}
		objNum++;
	}
	//Only time frame should be changed/modified
	frame_++;
}

void BrenderManager::add(shared_ptr<Brenderable> brenderable)
{
	brenderables_.push_back(brenderable);
}

void BrenderManager::setExportDir(string export_dir)
{
	EXPORT_DIR_ = export_dir;
}

