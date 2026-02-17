using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.CompilerServices;
using System.Runtime.InteropServices;
using System.Text;
using System.Threading.Tasks;

namespace DiffUrSolver
{
    public enum Methods
    {
        Euler,
        RK2,
        RK4,
        BackEuler
    }

    [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
    public delegate double F1(double x, double y);

    [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
    public delegate double F1System(double x, IntPtr y);

    [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
    public delegate double F(double x);

    [StructLayout(LayoutKind.Sequential)]
    public struct Task
    {
        public double x0;
        public double x1;
        public double y0;
        public double h;
        public IntPtr f;
        public Methods method;
    }

    [StructLayout(LayoutKind.Sequential)]
    public struct SecondTask
    {
        public double t0;
        public double t1;
        public double y_0;
        public double y1_0;
        public double h;
        public double A, B, C, D;
        public Methods method;
    }

    [StructLayout(LayoutKind.Sequential)]
    public struct SystemTask
    {
        public double t0;
        public double t1;
        public IntPtr y0;
        public uint n;
        public double h;
        public IntPtr f;
        public Methods method;
    }

    public class DiffUrSolver
    {

        public static IntPtr F1ToIntPtr(F1 func)
        {
            // Создаём массив IntPtr для хранения указателей на делегаты
            IntPtr functionPointer = new IntPtr();

            // Фиксируем делегаты в памяти, чтобы GC не перемещал их
            GCHandle handles = new GCHandle();

            handles = GCHandle.Alloc(func);
            functionPointer = Marshal.GetFunctionPointerForDelegate(func);

            GCHandle functionHandle = GCHandle.Alloc(functionPointer, GCHandleType.Pinned);

            return functionHandle.AddrOfPinnedObject();
        }

        public static IntPtr F1SystemToIntPtr(F1System[] system)
        {
            // Создаём массив IntPtr для хранения указателей на делегаты
            IntPtr[] functionPointers = new IntPtr[system.Length];

            // Фиксируем делегаты в памяти, чтобы GC не перемещал их
            GCHandle[] handles = new GCHandle[system.Length];



            for (int i = 0; i < system.Length; i++)
            {
                handles[i] = GCHandle.Alloc(system[i]);
                functionPointers[i] = Marshal.GetFunctionPointerForDelegate(system[i]);
            }


            GCHandle functionsArrayHandle = GCHandle.Alloc(functionPointers, GCHandleType.Pinned);

            return functionsArrayHandle.AddrOfPinnedObject();

        }

        public static IntPtr DoubleArrayToIntPtr(double[] arr)
        {
            int size = Marshal.SizeOf(typeof(double)) * arr.Length;
            IntPtr y0Ptr = Marshal.AllocHGlobal(size);
            Marshal.Copy(arr, 0, y0Ptr, arr.Length);

            return y0Ptr;

        }



        const string dllPath = "..\\..\\..\\..\\x64\\Release\\DiffUrSolverDLL.dll";

        [DllImport(dllPath, CallingConvention = CallingConvention.Cdecl)]
        public static extern int SolveDiffUr( Task task, ref double res);

        [DllImport(dllPath, CallingConvention = CallingConvention.Cdecl)]
        public static extern int SolveDiffUrAutoH(Task task, ref double res);

        [DllImport(dllPath, CallingConvention = CallingConvention.Cdecl)]
        public static extern int SolveDiffUrAutoHArr(Task task, double eps, double[] resX, double[] resY, ref int resSize);

        [DllImport(dllPath, CallingConvention = CallingConvention.Cdecl)]
        public static extern int SolveDiffUrAutoMethod(Task task, double eps, ref double res);

        [DllImport(dllPath, CallingConvention = CallingConvention.Cdecl)]
        public static extern int SolveDiffUrAutoMethodArr(Task task, double eps, double[] res);

        [DllImport(dllPath, CallingConvention = CallingConvention.Cdecl)]
        public static extern int SolveDiffUrAutoMethodArr2(Task task, double eps, double[] res, Methods[] methodsArr);

        [DllImport(dllPath, CallingConvention = CallingConvention.Cdecl)]
        public static extern int SolveDiffUrArr(Task task, double[] res);

        [DllImport(dllPath, CallingConvention = CallingConvention.Cdecl)]
        public static extern int SolveDiffUrSystem(SystemTask task, double[] res);

        [DllImport(dllPath, CallingConvention = CallingConvention.Cdecl)]
        public static extern int SolveDiffUrSystemArr(SystemTask task, double[] res);
    
        public static int SolveSecondDiffUr(SecondTask task,ref double res)
        {
            SystemTask systemTask = new SystemTask();

            systemTask.t0 = task.t0;
            systemTask.t1 = task.t1;
            systemTask.h = task.h;
            systemTask.method = task.method;
            systemTask.n = 2;

            double[] y0 = new double[6];
            y0[0] = task.y_0;
            y0[1] = task.y1_0;

            IntPtr y0Ptr = DoubleArrayToIntPtr(y0);
            systemTask.y0 = y0Ptr;


            double[] resultY = new double[systemTask.n];

            // Создаём массив IntPtr для хранения указателей на делегаты
            IntPtr[] functionPointers = new IntPtr[systemTask.n];

            // Фиксируем делегаты в памяти, чтобы GC не перемещал их
            GCHandle[] handles = new GCHandle[systemTask.n];

            F1System[] system = new F1System[systemTask.n];

            system[0] = delegate(double x, IntPtr y){
                double[] array = new double[2];
                Marshal.Copy(y, array, 0, 2);
                return array[1];
            };
            system[1] = delegate (double x, IntPtr y) {
                double[] array = new double[2];
                Marshal.Copy(y, array, 0, 2);
                return (task.D - task.B * array[1] - task.C * array[0]) / task.A;
            };

            for (int i = 0; i < system.Length; i++)
            {
                handles[i] = GCHandle.Alloc(system[i]);
                functionPointers[i] = Marshal.GetFunctionPointerForDelegate(system[i]);
            }


            GCHandle functionsArrayHandle = GCHandle.Alloc(functionPointers, GCHandleType.Pinned);

            systemTask.f = functionsArrayHandle.AddrOfPinnedObject();


            SolveDiffUrSystem(systemTask, resultY);
            res = resultY[0];
            return 0;
        }

        public static int SolveSecondDiffUrArr(SecondTask task, ref double[] res)
        {
            SystemTask systemTask = new SystemTask();

            systemTask.t0 = task.t0;
            systemTask.t1 = task.t1;
            systemTask.h = task.h;
            systemTask.method = task.method;
            systemTask.n = 2;

            double[] y0 = new double[6];
            y0[0] = task.y_0;
            y0[1] = task.y1_0;

            IntPtr y0Ptr = DoubleArrayToIntPtr(y0);
            systemTask.y0 = y0Ptr;

            F1System[] system = new F1System[systemTask.n];

            system[0] = delegate (double x, IntPtr y) {
                double[] array = new double[2];
                Marshal.Copy(y, array, 0, 2);
                return array[1];
            };

            system[1] = delegate (double x, IntPtr y) {
                double[] array = new double[2];
                Marshal.Copy(y, array, 0, 2);
                return (task.D - task.B * array[1] - task.C * array[0]) / task.A;
            };

            systemTask.f = F1SystemToIntPtr(system);

            double[] resultY = new double[10000];

            SolveDiffUrSystemArr(systemTask, resultY);

            double start = task.t0;
            int j = 0;

            res[0] = resultY[0];
            res[1] = resultY[1];

            while (task.t0 + j * task.h < task.t1)
            {
                j++;
                res[j * systemTask.n] =resultY[j*systemTask.n];
                res[j * systemTask.n + 1] = resultY[j * systemTask.n + 1];

                start = task.t0 + j* task.h;
            }

            j++;
            res[j * systemTask.n] = resultY[j * systemTask.n];
            res[j * systemTask.n + 1] = resultY[j * systemTask.n + 1];

            return 0;
        }
    }
}
