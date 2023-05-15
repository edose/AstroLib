using System.Runtime.InteropServices;
using AstroLib.Core;

// ReSharper disable IdentifierTypo
// ReSharper disable InconsistentNaming

namespace AstroLib.Core.SOFA;

public static partial class Sofa {
    
    /// <summary> Star-independent astrometry parameters, as a C++ struct, used to Marshal data to and
    /// to return data from SOFA C++ routines. Avoids multidimensional arrays, which can't be marshalled
    /// to and from unmanaged (C++) code. Data from C++ that are destined for .NET/C# multidimensional
    /// arrays must be recast from one-dimension to multidimensional after return from C++.
    ///  
    /// Vectors eb, eh, em, and v: with respect to axes of the BCRS (barycentric celestial ref. system).
    /// As defined in SOFA's sofa.h source file.</summary>
    [StructLayout(LayoutKind.Sequential)]
    public struct iauASTROM {
        public double pmt; // Proper motion time interval (Solar system barycenter, Julian years)

        [MarshalAs(UnmanagedType.ByValArray, SizeConst = 3)] // C++: double[3] eb
        public double[] eb; // Solar system barycenter to observer (vector, AU)

        [MarshalAs(UnmanagedType.ByValArray, SizeConst = 3)] // C++: double[3] eh
        public double[] eh; // Sun to observer (unit vector)

        public double em; // Distance from sun to observer (AU)

        [MarshalAs(UnmanagedType.ByValArray, SizeConst = 3)] // C++: double[3] v
        public double[] v; // Barycentric observer velocity (vector, units of c)

        public double bm1; // Reciprocal of Lorenz factor, i.e., sqrt(1-|v|^2)

        // [MarshalAs(UnmanagedType.ByValArray, SizeConst = 9)] // C++: double[3,3] bpn 2-dim as in SOFA
        // public double[,] bpn; // Bias-precession-nutation matrix
        
        [MarshalAs(UnmanagedType.ByValArray, SizeConst = 9)] // C++: double[9] bpn 1-dimensional
        public double[] bpn; // Bias-precession-nutation matrix

        public double along; // Longitude + s' + dETA(DUT) (radians)
        public double xp1; // Polar motion xp, with respect to local meridian (radians)
        public double yp1; // Polar motion yp, with respect to local meridian (radians)
        public double sphi; // Sine of geodetic latitude
        public double cphi; // Cosine of geodetic latitude
        public double diurab; // Magnitude of diurnal aberration vector
        public double era1; // "local" Earth rotation angle (radians)
        public double refa; // Refraction constant A (radians)
        public double refb; // Refraction constant B (radians)
    }

    /// <summary> Star-independent astrometry parameters, as a C#/.NET struct.
    /// This is the form exposed to C#, including multidimensional arrays recast from iauASTROM struct's
    /// one-dimensional (C++) arrays.
    /// 
    /// Vectors eb, eh, em, and v: with respect to axes of the BCRS (barycentric celestial ref. system).
    /// As defined in SOFA's sofa.h source file.</summary>
    [StructLayout(LayoutKind.Sequential)]
    public struct csharpASTROM {
        public double pmt; // Proper motion time interval (Solar system barycenter, Julian years)
        public double[] eb; // Solar system barycenter to observer (vector, AU)
        public double[] eh; // Sun to observer (unit vector)
        public double em; // Distance from sun to observer (AU)
        public double[] v; // Barycentric observer velocity (vector, units of c)
        public double bm1; // Reciprocal of Lorenz factor, i.e., sqrt(1-|v|^2)
        public double[,] bpn; // Bias-precession-nutation matrix  *** THIS IS RECAST FROM C++ one-dim.
        public double along; // Longitude + s' + dETA(DUT) (radians)
        public double xp1; // Polar motion xp, with respect to local meridian (radians)
        public double yp1; // Polar motion yp, with respect to local meridian (radians)
        public double sphi; // Sine of geodetic latitude
        public double cphi; // Cosine of geodetic latitude
        public double diurab; // Magnitude of diurnal aberration vector
        public double era1; // "local" Earth rotation angle (radians)
        public double refa; // Refraction constant A (radians)
        public double refb; // Refraction constant B (radians)
    }

    public static csharpASTROM TranslateToCsharp(iauASTROM cpp) {
        var csharp = new csharpASTROM();
        csharp.pmt = cpp.pmt;
        csharp.eb = new double[3];
        Array.Copy(cpp.eb, csharp.eb, cpp.eb.Length);
        csharp.eh = new double[3];
        Array.Copy(cpp.eh, csharp.eh, cpp.eh.Length);
        csharp.em = cpp.em;
        csharp.v = new double[3];
        Array.Copy(cpp.v, csharp.v, cpp.v.Length);
        csharp.bm1 = cpp.bm1;
        csharp.bpn = AstroMath.Reshape1dTo2dArray(cpp.bpn, 3, 3);
        csharp.along = cpp.along;
        csharp.xp1 = cpp.xp1;
        csharp.yp1 = cpp.yp1;
        csharp.sphi = cpp.sphi;
        csharp.cphi = cpp.cphi;
        csharp.diurab = cpp.diurab;
        csharp.era1 = cpp.era1;
        csharp.refa = cpp.refa;
        csharp.refb = cpp.refb;
        return csharp;
    }
    
    // /// <summary> Body parameters for light deflection.
    // /// As defined in SOFA's sofa.h source file.</summary>
    // [StructLayout(LayoutKind.Sequential)]
    // public struct iauLDBODY {
    //     public double bm; // Mass of the light-deflecting body (solar masses)
    //     public double dl; // Deflection limiter (radians^2/2)
    //
    //     [MarshalAs(UnmanagedType.ByValArray, SizeConst = 6)] //  C++: double[2,3] pv
    //     public double[,] pv; // Barycentric position/velocity of the light-deflecting body (AU, AU/day) 
    // }



}



