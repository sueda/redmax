#pragma once

#include <json/writer.h>
#include <json/value.h>
#include <json/json.h>

struct JSONwrapper
{
	Json::Value jsonscene;
	Json::Value frames;

	std::vector<Json::Value> blocks;
	std::vector<Json::Value> block_scales;
	std::vector<Json::Value> block_locations;
	std::vector<Json::Value> block_quats;

	JSONwrapper()
	{
		frames = Json::Value(Json::arrayValue);
	};
};