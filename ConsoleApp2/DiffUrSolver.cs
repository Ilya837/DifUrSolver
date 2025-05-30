using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.InteropServices;
using System.Text;
using System.Threading.Tasks;

namespace DiffUrSolver
{
    internal class DiffUrSolver
    {
        const string dllPath = "..\\..\\..\\..\\x64\\Release\\DiffUrSolverDLL.dll";

        [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
        public delegate double F1(double x, double y);

        [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
        public delegate double F1System(double x, double[] y);

        [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
        public delegate double StepF(double x0, double y0, double h, F1 f);

        [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
        public delegate double F(double x);

        [DllImport(dllPath, CallingConvention = CallingConvention.Cdecl)]
        public static extern int SolveDiffUr(
            double x0,
            double y0,
            double x1,
            double h,
            F1 f,
            ref double res,
            Methods method);

        [DllImport(dllPath, CallingConvention = CallingConvention.Cdecl)]
        public static extern int SolveDiffUrAutoH(double x0, double y0, double x1, double h, double eps,
            F1 f, ref double res, Methods method);

        [DllImport(dllPath, CallingConvention = CallingConvention.Cdecl)]
        public static extern int SolveDiffUrAutoHArr(double x0, double y0, double x1, double h, double eps,
            F1 f, double[] resX, double[] resY, ref int resSize, Methods method);

        [DllImport(dllPath, CallingConvention = CallingConvention.Cdecl)]
        public static extern int SolveDiffUrAutoMethod(double x0, double y0, double x1, double h, double eps,
            F1 f, ref double res, Methods method);

        [DllImport(dllPath, CallingConvention = CallingConvention.Cdecl)]
        public static extern int SolveDiffUrAutoMethodArr(double x0, double y0, double x1, double h, double eps,
            F1 f, double[] res, ref int sizeRes, Methods method);

        [DllImport(dllPath, CallingConvention = CallingConvention.Cdecl)]
        public static extern int SolveDiffUrAutoMethodArr2(double x0, double y0, double x1, double h, double eps,
            F1 f, double[] res, ref int sizeRes, Methods[] methodsArr, Methods method);

        [DllImport(dllPath, CallingConvention = CallingConvention.Cdecl)]
        public static extern int SolveDiffUrArr(double x0, double y0, double x1, double h, F1 f,
            double[] res, ref int sizeRes, Methods method);

        [DllImport(dllPath, CallingConvention = CallingConvention.Cdecl)]
        public static extern int SolveDiffUrSystem(double x0, double[] y0, double x1, uint n, double h, F1System[] f,
            double[] res, Methods method);

        [DllImport(dllPath, CallingConvention = CallingConvention.Cdecl)]
        public static extern int SolveDiffUrSystemArr(double x0, double[] y0, double x1, uint n, double h, F1System[] f,
            double[][] res, Methods method);
    }
}
