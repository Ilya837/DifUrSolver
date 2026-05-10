using System;
using System.Runtime.InteropServices;
using DiffUrSolver;
using OxyPlot;
using OxyPlot.Series;
using OxyPlot.Axes;
using OxyPlot.Core.Drawing;
using OxyPlot.WindowsForms;
using System.Threading.Tasks;
using System.Diagnostics;




class Program
{
    

    // Пример функции для решения
    public static double func1(double x, double y)
    {
        return Math.Cos(x);
    }

    public static double realFunc1(double x)
    {
        return Math.Sin(x);
    }

    public static double SampleFunction2(double x, double[] y)
    {
        return x * x + y[0];
    }

    public static double SampleFunction3(double x0, double y0, double h, DiffUrSolver.F1 f)
    {
        return 0;
    }

    public static double SampleFunction4(double x)
    {
        return x * x;
    }


    static void experiment1(DiffUrSolver.Methods method)
    {

        DiffUrSolver.Task task = new DiffUrSolver.Task();

        task.x0 = 0.0;
        task.y0 = 0.0;
        task.x1 = 1.0;
        task.h = 0.1;
        double result = 0;

        task.f = DiffUrSolver.DiffUrSolver.F1ToIntPtr(func1);
        task.method = method;

        int status = DiffUrSolver.DiffUrSolver.SolveDiffUr(task, ref result);

        Console.WriteLine(result.ToString());
        Console.WriteLine(realFunc1(task.x1));
    }

    static void experiment2(DiffUrSolver.Methods method)
    {

        DiffUrSolver.Task task = new DiffUrSolver.Task();

        task.x0 = 0.0;
        task.y0 = 0.0;
        task.x1 = 6.0;
        task.h = 0.3;
        double eps = 0.0001;
        double[] resultX = new double[1000];
        double[] resultY = new double[1000];
        int resultSize = 0;

        // Решение ОДУ методом Эйлера
        task.f = DiffUrSolver.DiffUrSolver.F1ToIntPtr(func1);
        task.method = method;

        int status = DiffUrSolver.DiffUrSolver.SolveDiffUrAutoHArr(task, eps, resultX, resultY, ref resultSize);

        //List<double> finalResult = result.TakeWhile(x => x != 0).ToList();




        var pm = new PlotModel();
        var ls = new LineSeries() { Title = "1", MarkerStroke = OxyColors.Black };
        var ls2 = new LineSeries() { Title = "2", MarkerStroke = OxyColors.Black };

        for (int i = 0; i < resultSize; i++)
        {
            //Console.WriteLine("{0:F7}|{1:F7}|{2:F7}|{3:F7}", resultX[i], resultY[i], RealSampleFunction(resultX[i]), RealSampleFunction(resultX[i]) - resultY[i]);

            ls.Points.Add(new DataPoint(resultX[i], resultY[i]));
            ls2.Points.Add(new DataPoint(resultX[i], realFunc1(resultX[i])));


        }

        pm.Series.Add(ls);
        pm.Series.Add(ls2);

    }
    static void TestMethodsEps(DiffUrSolver.F1 func, DiffUrSolver.F realFunc, DiffUrSolver.Methods method)
    {

        DiffUrSolver.Task task = new DiffUrSolver.Task();

        task.x0 = 0.0;
        task.y0 = 0.0;
        task.x1 = 4.0;
        task.h = 0.1;

        task.f = DiffUrSolver.DiffUrSolver.F1ToIntPtr(func);
        task.method = method;

        double[] result = new double[10000];
        double[] realResult = new double[10000];
        double localRes = 0;

        for (int i = 0; task.x0 + i * task.h <= task.x1; i++)
        {
            realResult[i] = realFunc(task.x0 + i * task.h);
        }

        DiffUrSolver.DiffUrSolver.SolveDiffUrArr(task, result);

        //for (int i = 0; task.x0 + i * task.h < task.x1; i++)
        //{
        //    Console.Write(result[i] + " ");
        //}

        //Console.WriteLine();

        //for (int i = 0; task.x0 + i * task.h < task.x1; i++)
        //{
        //    Console.Write(realResult[i] + " ");
        //}

        //Console.WriteLine();

        double globalError1 = -10000000;
        double localError1 = -10000000;

        double startx0 = task.x0;
        double starty0 = task.y0;
        double startx1 = task.x1;

        for (int i = 0; task.x0 + task.h < startx1; i++)
        {

            task.x0 = startx0 + task.h * i;
            task.y0 = realFunc(task.x0);
            task.x1 = task.x0 + task.h;

            if (Math.Abs(result[i] - realResult[i]) > globalError1)
            {
                globalError1 = Math.Abs(result[i] - realResult[i]);
            }

            DiffUrSolver.DiffUrSolver.SolveDiffUr(task, ref localRes);

            if (Math.Abs(localRes - realResult[i + 1]) > localError1)
            {
                localError1 = Math.Abs(localRes - realResult[i + 1]);
            }
        }

        task.x0 = startx0;
        task.y0 = starty0;
        task.x1 = startx1;
        task.h /= 2;

        for (int i = 0; task.x0 + i * task.h <= task.x1; i++)
        {
            realResult[i] = realFunc(task.x0 + i * task.h);
        }

        DiffUrSolver.DiffUrSolver.SolveDiffUrArr(task, result);

        //for (int i = 0; task.x0 + i * task.h < task.x1; i++)
        //{
        //    Console.Write(result[i] + " ");
        //}

        //Console.WriteLine();

        //for (int i = 0; task.x0 + i * task.h < task.x1; i++)
        //{
        //    Console.Write(realResult[i] + " ");
        //}

        //Console.WriteLine();

        double globalError2 = -10000000;
        double localError2 = -10000000;

        startx0 = task.x0;
        startx1 = task.x1;
        starty0 = task.y0;

        for (int i = 0; task.x0 + task.h < startx1; i++)
        {

            task.x0 = startx0 + task.h * i;
            task.y0 = realFunc(task.x0);
            task.x1 = task.x0 + task.h;

            if (Math.Abs(result[i] - realResult[i]) > globalError2)
            {
                globalError2 = Math.Abs(result[i] - realResult[i]);
            }

            DiffUrSolver.DiffUrSolver.SolveDiffUr(task, ref localRes);

            if (Math.Abs(localRes - realResult[i + 1]) > localError2)
            {
                localError2 = Math.Abs(localRes - realResult[i + 1]);
            }
        }

        task.x0 = startx0;
        task.y0 = starty0;
        task.x1 = startx1;
        task.h *= 2;

        Console.WriteLine("          |    Global error     |   Local error");

        Console.WriteLine(" h = {0:F2} |     {1:F10}    |   {2:F10}", task.h, globalError1, localError1);

        Console.WriteLine(" h = {0:F2} |     {1:F10}    |   {2:F10}", task.h / 2, globalError2, localError2);

        Console.WriteLine("          |     h^{0:F10}  | h^{1:F10}", Math.Log2(globalError1 / globalError2), Math.Log2(localError1 / localError2));


    }

    static void FullTestMethodsEps(DiffUrSolver.F1 func, DiffUrSolver.F realFunc)
    {
        Console.WriteLine("                      Euler");
        TestMethodsEps(func, realFunc, DiffUrSolver.Methods.Euler);
        Console.WriteLine();
        Console.WriteLine("                       RK2");
        TestMethodsEps(func, realFunc, DiffUrSolver.Methods.RK2);
        Console.WriteLine();
        Console.WriteLine("                       RK4");
        TestMethodsEps(func, realFunc, DiffUrSolver.Methods.RK4);
        Console.WriteLine();
    }

    public static double SystemFunctionSingle(double x, IntPtr y)
    {

        double[] array = new double[1];
        Marshal.Copy(y, array, 0, 1);

        return Math.Cos(array[0]);
    }
    public static double SystemFunction1(double x, IntPtr y)
    {

        double[] array = new double[2];
        Marshal.Copy(y, array, 0, 2);


        return array[0] - 5 * array[1];
    }

    public static double SystemFunction2(double x, IntPtr y)
    {
        double[] array = new double[2];
        Marshal.Copy(y, array, 0, 2);

        return 5 * array[0] + array[1];
    }

    public static double RealSystemFunction1(double x)
    {
        return Math.Pow(Math.E, x) * (Math.Sin(-5 * x) + Math.Cos(-5 * x));
    }

    public static double RealSystemFunction2(double x)
    {
        return Math.Pow(Math.E, x) * (Math.Sin(5 * x) + Math.Cos(5 * x));
    }

    static void experiment3(DiffUrSolver.Methods method)
    {
        DiffUrSolver.SystemTask task = new DiffUrSolver.SystemTask();

        task.t0 = 0;
        task.t1 = 4.0;
        task.h = 0.001;
        task.method = method;
        task.n = 2;

        double[] y0 = new double[6];
        y0[0] = 1;
        y0[1] = 1;

        int size = Marshal.SizeOf(typeof(double)) * y0.Length;
        IntPtr y0Ptr = Marshal.AllocHGlobal(size);
        Marshal.Copy(y0, 0, y0Ptr, y0.Length);

        task.y0 = y0Ptr;


        double[] resultY = new double[task.n];

        // Создаём массив IntPtr для хранения указателей на делегаты
        IntPtr[] functionPointers = new IntPtr[2];

        // Фиксируем делегаты в памяти, чтобы GC не перемещал их
        GCHandle[] handles = new GCHandle[2];

        F1System[] system = new F1System[2];

        system[0] = SystemFunction1;
        system[1] = SystemFunction2;

        for (int i = 0; i < system.Length; i++)
        {
            handles[i] = GCHandle.Alloc(system[i]);
            functionPointers[i] = Marshal.GetFunctionPointerForDelegate(system[i]);
        }


        GCHandle functionsArrayHandle = GCHandle.Alloc(functionPointers, GCHandleType.Pinned);

        task.f = functionsArrayHandle.AddrOfPinnedObject();


        DiffUrSolver.DiffUrSolver.SolveDiffUrSystem(task, resultY);

        for (int i = 0; i < task.n; i++)
        {
            Console.WriteLine(resultY[i]);
        }

        Console.WriteLine(RealSystemFunction1(task.t1));
        Console.WriteLine(RealSystemFunction2(task.t1));

    }
    static void experiment4(DiffUrSolver.Methods method)
    {
        DiffUrSolver.SystemTask task = new DiffUrSolver.SystemTask();

        task.t0 = 0;
        task.t1 = 4.0;
        task.h = 0.001;
        task.method = method;
        task.n = 1;

        double[] y0 = new double[1];
        y0[0] = 0;
        //y0[1] = 1;

        task.y0 = DiffUrSolver.DiffUrSolver.DoubleArrayToIntPtr(y0);

        double[] resultY = new double[task.n];

        F1System[] system = new F1System[task.n];

        system[0] = SystemFunctionSingle;

        task.f = DiffUrSolver.DiffUrSolver.F1SystemToIntPtr(system);

        DiffUrSolver.DiffUrSolver.SolveDiffUrSystem(task, resultY);

        for (int i = 0; i < task.n; i++)
        {
            Console.WriteLine(resultY[i]);
        }

        Console.WriteLine(realFunc1(task.t1));

    }

    static void experiment5(DiffUrSolver.Methods method)
    {
        double c = 1.0;
        double m = 1.0;
        double k = 1.0;

        DiffUrSolver.SecondTask task = new DiffUrSolver.SecondTask();
        task.t0 = 0;
        task.t1 = 10.0;
        task.y_0 = 2;
        task.y1_0 = 0;
        task.h = 0.1;
        task.A = 1;
        task.B = c / m;
        task.C = k / m;
        task.D = 0;
        task.method = method;
        double result = 0;
        DiffUrSolver.DiffUrSolver.SolveSecondDiffUr(task, ref result);
        Console.WriteLine(result.ToString());
    }
    static void experiment5Arr(DiffUrSolver.Methods method)
    {
        double c = 4.47;
        double m = 10.0;
        double k = 50.0;

        DiffUrSolver.SecondTask task = new DiffUrSolver.SecondTask();
        task.t0 = 0;
        task.t1 = 10;
        task.y_0 = 0.01;
        task.y1_0 = 0;
        task.h = 0.1;
        task.A = 1;
        task.B = c / m;
        task.C = k / m;
        task.D = 0;
        task.method = method;
        double[] result = new double[10000];
        DiffUrSolver.DiffUrSolver.SolveSecondDiffUrArr(task, ref result);


        double start = task.t0;
        int j = 0;

        while (start < task.t1)
        {
            Console.WriteLine(result[j * 2]);
            j++;
            start = task.t0 + j * task.h;
        }
    }

    struct BrusselatorSystem
    {
        public IntPtr f;
        public IntPtr y0;
        public double t0;
        public int status;
    }

    static BrusselatorSystem CreateBrusselatorSystem(int n, double a)
    {

        BrusselatorSystem res = new BrusselatorSystem();

        if (n < 1)
        {

            res.status = -1;
            return res;
        }

        F1System[] system = new F1System[n * 2];
        double[] y = new double[2 * n];

        for (int i = 0; i < n; i++)
        {

            int i1 = (n + i - 1) % n;
            int i2 = i;
            int i3 = (n + i + 1) % n;

            F1System u = (double x, IntPtr y) =>
            {
                double[] array = new double[2 * n];
                Marshal.Copy(y, array, 0, 2 * n);

                return 1 + array[i2 * 2] * array[i2 * 2] * array[i2 * 2 + 1] -
                4 * array[i2 * 2] + a * (n + 1) * (n + 1) *
                (array[i1 * 2] - 2 * array[i2 * 2] + array[i3 * 2]);
            };

            F1System v = (double x, IntPtr y) =>
            {
                double[] array = new double[2 * n];
                Marshal.Copy(y, array, 0, 2 * n);

                return 3 * array[i2 * 2] - array[i2 * 2] * array[i2 * 2] * array[i2 * 2 + 1] +
                a * (n + 1) * (n + 1) *
                (array[i1 * 2 + 1] - 2 * array[i2 * 2 + 1] + array[i3 * 2 + 1]);
            };

            system[i * 2] = u;
            system[i * 2 + 1] = v;

            y[i * 2] = 1 + Math.Sin(2 * Math.PI * i / n);
            y[i * 2 + 1] = 3;

        }

        res.status = 0;
        res.f = DiffUrSolver.DiffUrSolver.F1SystemToIntPtr(system);
        res.y0 = DiffUrSolver.DiffUrSolver.DoubleArrayToIntPtr(y);
        res.t0 = 0;


        return res;

    }

    static void experiment6(double h, Methods method,uint treadCount)
    {
        DiffUrSolver.DiffUrSolver.SetThreadCount(treadCount);
        DiffUrSolver.SystemTask systemTask = new DiffUrSolver.SystemTask();

        uint n = 500;

        BrusselatorSystem b = CreateBrusselatorSystem(((int)n), 0.02);


        systemTask.y0 = b.y0;
        systemTask.f = b.f;
        systemTask.t0 = b.t0;
        systemTask.t1 = 10;
        systemTask.h = h;
        systemTask.n = 2 * n;
        systemTask.method = method;

        int size1 = (int)Math.Round((systemTask.t1 - systemTask.t0) / systemTask.h) + 10;
        //double[] resultX = new double[size1];
        double[] resultY = new double[systemTask.n];

        Stopwatch stopwatch = new Stopwatch();

        DiffUrSolver.DiffUrSolver.InitMemory(systemTask.n);

        stopwatch.Start();

        int status = DiffUrSolver.DiffUrSolver.SolveDiffUrSystem(systemTask, resultY);

        stopwatch.Stop();

        System.Console.WriteLine(stopwatch.ElapsedMilliseconds);
        System.Console.WriteLine("result:");

        for (int i = 0; i < systemTask.n; i++)
        {
            Console.WriteLine(resultY[i]);
        }

    }




    static void Main(string[] args)
    {
        //experiment1(Methods.Euler);
        //experiment1(Methods.BackEuler);

        //experiment5Arr(Methods.RK4);

        //experiment2(DiffUrSolver.Methods.RK4);

        //FullTestMethodsEps(func1,realFunc1);

        //FullTestMethodsEps(func1,realFunc1);

        //experiment4(Methods.Euler);

        experiment6(0.00001, Methods.Euler, 1);

        //experiment6(0.02, Methods.BackEuler, 4);

        //for (uint i = 1; i <= 4; i++) {
        //    experiment6(0.02, Methods.BackEuler, i);
        //}
    }
}
