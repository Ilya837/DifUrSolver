using System.Windows.Forms;
using DiffUrSolver;
using OxyPlot;
using OxyPlot.Annotations;
using OxyPlot.Series;
using OxyPlot.WindowsForms;

namespace WinFormsApp
{
    public partial class Form2 : Form
    {
        public Form2()
        {
            InitializeComponent();
            DrawAutoHPlot(SampleFunction, RealSampleFunction,DiffUrSolver.Methods.RK4);
            //DrawAutoMethodPlot(SampleFunction, RealSampleFunction, DiffUrSolver.Methods.Euler);
        }

        public static double SampleFunction(double x, double y)
        {
            return Math.Cos(x);
        }

        public static double RealSampleFunction(double x)
        {
            return Math.Sin(x);
        }

        private const AnnotationLayer HiddenLayer = (AnnotationLayer)100;
        private const AnnotationLayer VisibleLayer = AnnotationLayer.AboveSeries;

        private class EpsAnnotation : TextAnnotation { };

        private class MethodAnnotation : TextAnnotation { };

        private void DrawAutoHPlot(DiffUrSolver.F1 func, DiffUrSolver.F realFunc,DiffUrSolver.Methods method)
        {
            var plotModel = new PlotModel { Title = "Пример простого графика" };

            var lineSeries = new LineSeries
            {
                Title = "My",
                MarkerType = MarkerType.Circle,
                MarkerSize = 3,
                MarkerStroke = OxyColors.Red,
                MarkerFill = OxyColors.Red,
                
            };

            var lineSeries2 = new LineSeries
            {
                Title = "Real",
                MarkerType = MarkerType.Circle,
                MarkerSize = 3,
                MarkerStroke = OxyColors.Green,
                MarkerFill= OxyColors.Green,
            };

            double x0 = 0.0;
            double y0 = 0.0;
            double x1 = 6.0;
            double h = 0.3;
            double eps = 0.0001;
            double[] resultX = new double[1000];
            double[] resultY = new double[1000];
            int resultSize = 0;



            int status = DiffUrSolver.DiffUrSolver.SolveDiffUrAutoHArr(x0, y0, x1, h, eps, func, resultX, resultY, ref resultSize, method);


            for(int i =0; i< resultSize; i++)
            {
                lineSeries.Points.Add(new DataPoint(resultX[i], resultY[i]));
                lineSeries2.Points.Add(new DataPoint(resultX[i], realFunc(resultX[i])));
                double diff = Math.Abs(resultY[i] - realFunc(resultX[i]));
                plotModel.Annotations.Add(new EpsAnnotation
                {
                    Text = diff.ToString("F6"),
                    TextPosition = new DataPoint(resultX[i], (resultY[i] + realFunc(resultX[i])) / 2),
                    Background = OxyColors.White,
                    TextColor = OxyColors.Black,
                    Layer = HiddenLayer

                });
            }


            plotModel.Series.Add(lineSeries);
            plotModel.Series.Add(lineSeries2);

            var plotView = new PlotView
            {
                Dock = DockStyle.Fill,
                Model = plotModel
            };

            this.Controls.Add(plotView);

            this.KeyPreview = true;
            this.KeyDown += (sender, e) =>
            {
                if (e.Control && e.KeyCode == Keys.E)
                {
                    foreach (var annotation in plotModel.Annotations.OfType<EpsAnnotation>())
                    {
                        annotation.Layer = annotation.Layer == HiddenLayer ? VisibleLayer : HiddenLayer;
                    }
                    plotModel.InvalidatePlot(true);
                }

                if (e.Control && e.KeyCode == Keys.M)
                {
                    foreach (var annotation in plotModel.Annotations.OfType<MethodAnnotation>())
                    {
                        annotation.Layer = annotation.Layer == HiddenLayer ? VisibleLayer : HiddenLayer;
                    }
                    plotModel.InvalidatePlot(true);
                }
            };

        }

        private void DrawAutoMethodPlot(DiffUrSolver.F1 func, DiffUrSolver.F realFunc, DiffUrSolver.Methods method)
        {
            var plotModel = new PlotModel { Title = "Пример простого графика" };

            var lineSeries = new LineSeries
            {
                Title = "My",
                MarkerType = MarkerType.Circle,
                MarkerSize = 3,
                MarkerStroke = OxyColors.Red,
                MarkerFill = OxyColors.Red,

            };

            var lineSeries2 = new LineSeries
            {
                Title = "Real",
                MarkerType = MarkerType.Circle,
                MarkerSize = 3,
                MarkerStroke = OxyColors.Green,
                MarkerFill = OxyColors.Green,
            };

            double x0 = 0.0;
            double y0 = 0.0;
            double x1 = 6.0;
            double h = 0.05;
            double eps = 0.0001;
            double[] resultX = new double[1000];
            double[] resultY = new double[1000];
            Methods[] methodsArr = new Methods[1000];
            int resultSize = 0;



            int status = DiffUrSolver.DiffUrSolver.SolveDiffUrAutoMethodArr2(x0, y0, x1, h, eps, func, resultY, ref resultSize, methodsArr, method);


            lineSeries.Points.Add(new DataPoint(x0 , resultY[0]));
            lineSeries2.Points.Add(new DataPoint(x0, realFunc(x0)));
            double diff = Math.Abs(resultY[0] - realFunc(x0));
            plotModel.Annotations.Add(new EpsAnnotation
            {
                Text = diff.ToString("F6"),
                TextPosition = new DataPoint(x0, (resultY[0] + realFunc(x0)) / 2),
                Background = OxyColors.White,
                TextColor = OxyColors.Black,
                Layer = HiddenLayer

            });

            for (int i = 1; i < resultSize; i++)
            {
                lineSeries.Points.Add(new DataPoint(x0 + i * h, resultY[i]));
                lineSeries2.Points.Add(new DataPoint(x0 + i * h, realFunc(x0 + i * h)));
                diff = Math.Abs(resultY[i] - realFunc(x0 + i * h));
                plotModel.Annotations.Add(new EpsAnnotation
                {
                    Text = diff.ToString("F6"),
                    TextPosition = new DataPoint(x0 + i * h, (resultY[i] + realFunc(x0 + i * h)) / 2),
                    Background = OxyColors.White,
                    TextColor = OxyColors.Black,
                    Layer = HiddenLayer

                });

                

                plotModel.Annotations.Add(new MethodAnnotation
                {
                    Text = methodsArr[i].ToString(),
                    TextPosition = new DataPoint((x0 + i * h + x0 + (i-1) * h)/2, (resultY[i] + resultY[i-1]) / 2),
                    Background = OxyColors.White,
                    TextColor = OxyColors.Black,
                    Layer = HiddenLayer

                });
                
            }


            plotModel.Series.Add(lineSeries);
            plotModel.Series.Add(lineSeries2);

            var plotView = new PlotView
            {
                Dock = DockStyle.Fill,
                Model = plotModel
            };

            this.Controls.Add(plotView);

            this.KeyPreview = true;
            this.KeyDown += (sender, e) =>
            {
                if (e.Control && e.KeyCode == Keys.E)
                {
                    foreach (var annotation in plotModel.Annotations.OfType<EpsAnnotation>())
                    {
                        annotation.Layer = annotation.Layer == HiddenLayer ? VisibleLayer : HiddenLayer;
                    }
                    plotModel.InvalidatePlot(true);
                }

                if (e.Control && e.KeyCode == Keys.M)
                {
                    foreach (var annotation in plotModel.Annotations.OfType<MethodAnnotation>())
                    {
                        annotation.Layer = annotation.Layer == HiddenLayer ? VisibleLayer : HiddenLayer;
                    }
                    plotModel.InvalidatePlot(true);
                }
            };

        }
    }
}