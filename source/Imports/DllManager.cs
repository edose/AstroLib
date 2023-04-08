namespace AstroLib.Imports; 
using System.Runtime.InteropServices;

public static class DllManager {

    [DllImport("kernel32.dll", CharSet = CharSet.Unicode, SetLastError = true)]
    private static extern IntPtr LoadLibrary(string libraryName);

    [DllImport("kernel32.dll", CharSet = CharSet.Unicode, SetLastError = true)]
    [return: MarshalAs(UnmanagedType.Bool)]
    private static extern bool SetDllDirectory(string lpPathName);

    private static readonly object DllLockObject = new object();

    public static void LoadDllFromRelativePath(string dllRelativePath) {
        var baseDirectory = AppDomain.CurrentDomain.BaseDirectory;
        var dllFullpath = Path.Combine(baseDirectory, "External", "x64", dllRelativePath);
        LoadDllFromFullpath(dllFullpath);            
    }

    public static void LoadDllFromFullpath(string dllFullpath) {
        var dllDirectory = Path.GetDirectoryName(dllFullpath)!;
        lock (DllLockObject) {
            SetDllDirectory(dllDirectory);
            if (LoadLibrary(dllFullpath) == IntPtr.Zero) {
                var error = Marshal.GetLastWin32Error().ToString().Trim();
                var message = $"DllLoader failed to load library {dllFullpath}, error code >{error}<";
                Console.WriteLine($">>> {message}");
            }
            SetDllDirectory("");
        }
    }
}