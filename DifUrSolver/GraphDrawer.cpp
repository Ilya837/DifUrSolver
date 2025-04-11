#include "string"
#include "GraphDrawer.h"
#include "discpp.h"

void GraphDrawer::DrawGraph(double* x, double* y, ui* n,ui k,double xmin, double xmax, double xstep, double ymin, double ymax, double ystep,std::string* colors, std::string* legend)
{
    Dislin g;
    int ic;

    g.metafl("cons");
    g.scrmod("revers");
    g.disini();
    g.pagera();
    g.complx();
    g.axspos(300, 1800);
    g.axslen(2600, 1750);

    g.name("X-axis", "x");
    g.name("Y-axis", "y");

    g.labdig(-1, "x");
    g.ticks(9, "x");
    g.ticks(10, "y");

    /*g.titlin("Sravnenie", 1);
    g.titlin("Euler, RealF", 3);*/

    ic = g.intrgb(0.95, 0.95, 0.95);
    g.axsbgd(ic);

    g.graf(xmin, xmax, xmin, xstep, ymin, ymax, ymin, ystep);
    g.setrgb(0.7, 0.7, 0.7);
    g.grid(1, 1);

    g.color("fore");
    g.height(50);
    g.title();

    int sum = 0;
    for (ui i = 0; i < k; i++) {
        g.color(colors[i].c_str());
            g.curve(x + sum, y + sum, n[i]);
        sum += n[i];
    }
    g.disfin();
}
