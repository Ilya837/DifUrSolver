#pragma once

#define ui unsigned int

class GraphDrawer
{
public:
	void static DrawGraph(double* x, double* y, ui* n, ui k, double xmin, double xmax, double xstep, double ymin, double ymax, double ystep, std::string* colors, std::string* legend);
};

