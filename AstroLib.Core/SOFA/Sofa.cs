using System.Runtime.InteropServices;
using AstroLib.Core.Utils;

// ReSharper disable StringLiteralTypo
// ReSharper disable CommentTypo
// ReSharper disable IdentifierTypo
// ReSharper disable InconsistentNaming

namespace AstroLib.Core.SOFA;

/*########################################################################################################*/
// Due diligence statement from the principal author of this class and its resulting assembly (.dll file):
// This work (AstroLib.Sofa):
// (i) calls routines and computations that I derived from software provided by SOFA under license; and
// (ii) does not itself constitute software provided by and/or endorsed by SOFA.
// (iii) calls *compiled* C-language functions from C-language source code provided by SOFA,
//       but does not distribute that original source code. We offer instead a .NET-based API to SOFA
//       functions, which API is especially suitable for calling from another user's own C#-language code.
/*########################################################################################################*/
//  Correspondence concerning SOFA software should be addressed as follows:
//
//      By email:  sofa@ukho.gov.uk
//      By post:   IAU SOFA Center
//                 HM Nautical Almanac Office
//                 UK Hydrographic Office
//                 Admiralty Way, Taunton
//                 Somerset, TA1 2DN
//                 United Kingdom
/*########################################################################################################*/

/// <summary>Class Sofa: Low-level .NET representations of selected IAU SOFA functions.</summary>
public static class Sofa {
    private const string DllDirectory = "C:/DevCS/SOFA_dll/x64/Release";
    private const string DllFilename  = "SOFA_dll.dll";

    static Sofa() {
        var dllFullpath = Path.Combine(DllDirectory, DllFilename);
        DllManager.LoadDllFromFullpath(dllFullpath);
    }

    /// <summary> Star-independent astrometry parameters, as a C++ struct.
    /// Vectors eb, eh, em, and v: with respect to axes of the BCRS (barycentric celestial ref. system).
    /// As defined in SOFA's sofa.h source file.</summary>
    [StructLayout(LayoutKind.Sequential)]
    public struct iauASTROM {
        public double pmt;     // Proper motion time interval (Solar system barycenter, Julian years)
        [MarshalAs(UnmanagedType.ByValArray, SizeConst = 3)]  // C++: double[3] eb
        public double[] eb;    // Solar system barycenter to observer (vector, AU)
        [MarshalAs(UnmanagedType.ByValArray, SizeConst = 3)]  // C++: double[3] eh
        public double[] eh;    // Sun to observer (unit vector)
        public double em;      // Distance from sun to observer (AU)
        [MarshalAs(UnmanagedType.ByValArray, SizeConst = 3)]  // C++: double[3] v
        public double[] v;     // Barycentric observer velocity (vector, units of c)
        public double bm1;     // Reciprocal of Lorenz factor, i.e., sqrt(1-|v|^2)
        [MarshalAs(UnmanagedType.ByValArray, SizeConst = 9)]  // C++: double[3,3] bpn
        public double[,] bpn;  // Bias-precession-nutation matrix
        public double along;   // Longitude + s' + dETA(DUT) (radians)
        public double xp1;     // Polar motion xp, with respect to local meridian (radians)
        public double yp1;     // Polar motion yp, with respect to local meridian (radians)
        public double sphi;    // Sine of geodetic latitude
        public double cphi;    // Cosine of geodetic latitude
        public double diurab;  // Magnitude of diurnal aberration vector
        public double era1;    // "local" Earth rotation angle (radians)
        public double refa;    // Refraction constant A (radians)
        public double refb;    // Refraction constant B (radians)
    }

    /// <summary> Body parameters for light deflection.
    /// As defined in SOFA's sofa.h source file.</summary>
    [StructLayout(LayoutKind.Sequential)]
    public struct iauLDBODY {
        public double bm;     // Mass of the light-deflecting body (solar masses)
        public double dl;     // Deflection limiter (radians^2/2)
        [MarshalAs(UnmanagedType.ByValArray, SizeConst = 6)]  //  C++: double[2,3] pv
        public double[,] pv;  // Barycentric position/velocity of the light-deflecting body (AU, AU/day) 
    }
    
    
#region "Astronomy/Calendars"
    /*######### Astronomy/Calendars ######################################################################*/
    // No need to expose as public *for now* (reproduce these as needed, in C# within AstroLib):
    //      iauEpb    Julian Date to Besselian Epoch.
    //      iauEpb2jd Besselian Epoch to Julian Date.
    //      iauEpj    Julian Date to Julian Epoch.
    //      iauEpj2jd Julian Epoch to Julian Date.
    //      iauJdcalf Julian Date to Gregorian Calendar, expressed in a form convenient
    //                    for formatting messages: rounded to a specified precision.
    // .....................................................................................................

    /// <summary>Sofa.Cal2jd(): Gregorian Calendar to Julian Date.</summary>
    /// <param name="iy">Year in Gregorian calendar</param>
    /// <param name="im">Month in Gregorian calendar</param>
    /// <param name="id">Day in Gregorian calendar</param>
    /// <returns>2-tuple of doubles:
    ///              djm0: MJD zero−point, always 2400000.5,
    ///              djm:  Modified Julian Date for 0 hrs. </returns>
    public static Tuple<double, double> Cal2jd(int iy, int im, int id) {
        double djm0 = 0, djm = 0;
        S_Cal2jd(iy, im, id, ref djm0, ref djm);
        return Tuple.Create(djm0, djm); }
    [DllImport(DllFilename, EntryPoint = "iauCal2jd", CallingConvention = CallingConvention.Cdecl)]
    static extern double S_Cal2jd(int iy, int im, int id, ref double djm0, ref double djm);
    
    /// <summary>Sofa.Jd2cal(): Julian Date to Gregorian year, month, day, and fraction of a day.</summary>
    /// <param name="dj1">Part 1 of 2-part Julian Date</param>
    /// <param name="dj2">Part 2 of 2-part Julian Date</param>
    /// <returns>4-tuple:
    ///              int iy, im, and id: year, month, and day,
    ///              double fd: fraction of day.</returns>
    public static Tuple<int, int, int, double> Jd2cal(double dj1, double dj2) {
        int iy = 0, im = 0, id = 0;
        double fd = 0;
        S_Jd2cal(dj1, dj2, ref iy, ref im, ref id, ref fd);
        return Tuple.Create(iy, im, id, fd); }
    [DllImport(DllFilename, EntryPoint = "iauJd2cal", CallingConvention = CallingConvention.Cdecl)]
    static extern double S_Jd2cal(double dj1, double dj2, 
        ref int iy, ref int im, ref int id, ref double fd);    
#endregion
    
#region "Astronomy/Astrometry"    
    /*######### Astronomy/Astrometry #####################################################################*/
    // No need to expose as public *for now* (reproduce these as needed, in C# within AstroLib):
    //      iauApcg   For a geocentric observer, prepare star−independent astrometry parameters
    //                    for transformations between ICRS and GCRS coordinates.
    //                    The Earth ephemeris is supplied by the caller.
    //      *** iauApcg13 As iauApcg, except caller supplies the date, and SOFA models used to predict
    //                    the Earth ephemeris.
    //      iauApci   For a terrestrial observer, prepare star−independent astrometry parameters
    //                    for transformations between ICRS and geocentric CIRS coordinates.
    //                    The Earth ephemeris and CIP/CIO are supplied by the caller.
    //      iauApci13 As iauApci, except caller supplies the date, and SOFA models used to predict
    //                    the Earth ephemeris and CIP/CIO.
    //      iauApco   For a terrestrial observer, prepare star−independent astrometry parameters
    //                    for transformations between ICRS and observed coordinates. The caller supplies
    //                    the Earth ephemeris, the Earth rotation information and the refraction
    //                    constants as well as the site coordinates.
    //      iauApco13 As iauApco, except caller supplies UTC, site coordinates, ambient air conditions
    //                    and observing wavelength, and SOFA models are to obtain
    //                    the Earth ephemeris, CIP/CIO and refraction constants.
    //      iauApcs   For an observer whose geocentric position and velocity are known, prepare
    //                    star−independent astrometry parameters for transformations between ICRS
    //                    and GCRS. The Earth ephemeris is supplied by the caller.
    //      iauApcs13 As iauApcs, except Earth ephemeris is from SOFA models.
    //      iauAper   In the star−independent astrometry parameters, update only the
    //                    Earth rotation angle, supplied by the caller explicitly.
    //      iauAper13 As iauAper, except caller provides UT1, (n.b. not UTC).
    //      iauApio   For a terrestrial observer, prepare star−independent astrometry parameters
    //                    for transformations between CIRS and observed coordinates. The caller
    //                    supplies the Earth orientation information and the refraction constants
    //                    as well as the site coordinates.
    //      iauApio13 As iauApio, except caller supplies UTC, site coordinates, ambient air conditions
    //                    and observing wavelength.
    //      iauAtccq  Quick transformation of a star’s ICRS catalog entry (epoch J2000.0)
    //                    into ICRS astrometric place, given precomputed star−independent
    //                    astrometry parameters.
    //      iauAtciq  Quick ICRS, epoch J2000.0, to CIRS transformation, given precomputed
    //                    star−independent astrometry parameters.
    //      iauAtciqn Quick ICRS, epoch J2000.0, to CIRS transformation, given precomputed
    //                    star−independent astrometry parameters plus a list of light−deflecting bodies.
    //      iauAtciqz Quick ICRS to CIRS transformation, given precomputed star−independent
    //                    astrometry parameters, and assuming zero parallax and proper motion.
    //      iauAticq  Quick CIRS RA,Dec to ICRS astrometric place, given the star−independent
    //                    astrometry parameters.
    //      iauAticqn Quick CIRS to ICRS astrometric place transformation, given the star−independent
    //                    astrometry parameters plus a list of light−deflecting bodies.
    //      iauAtioq  Quick CIRS to observed place transformation.
    //      iauAtoiq  Quick observed place to CIRS, given the star−independent astrometry parameters.
    //      iauLd     Apply light deflection by a solar−system body, as part of transforming
    //                    coordinate direction into natural direction.
    //      iauLdn    For a star, apply light deflection by multiple solar−system bodies, as part
    //                    of transforming coordinate direction into natural direction.
    //      iauLdsun  Deflection of starlight by the Sun.
    //      iauPmpx   Proper motion and parallax.
    // .....................................................................................................
    
    // .....................................................................................................
    // From SOFA C Manual, release 18:
    // ...several functions...insert...into the astrom [iauASTROM] structure
    // star−independent parameters needed for the chain of astrometric transformations
    // ICRS <−> GCRS <−> CIRS <−> observed.
    //
    // The various functions support different classes of observer and portions
    // of the transformation chain:
    //
    //    ____functions____     observer     transformation
    // 
    //   iauApcg   iauApcg13   geocentric    ICRS <−> GCRS
    //   iauApci   iauApci13   terrestrial   ICRS <−> CIRS
    //   iauApco   iauApco13   terrestrial   ICRS <−> observed
    //   iauApcs   iauApcs13   space         ICRS <−> GCRS
    //   iauAper   iauAper13   terrestrial   update Earth rotation
    //   iauApio   iauApio13   terrestrial   CIRS <−> observed
    // .....................................................................................................
    
    
    /// <summary> Sofa.Ab(): Apply aberration to transform natural direction into proper direction.
    /// </summary>
    /// <param name="pnat">Natural direction to the source (unit vector).</param>
    /// <param name="v">Observer barycentric velocity in units of c.</param>
    /// <param name="s">Distance between the Sun and the observer (AU).</param>
    /// <param name="bm1">Sqrt(1 - |v|^2), the reciprocal of Lorenz factor.</param>
    /// <returns>Proper direction to the source (unit vector).</returns>
    public static double[] Ab(double[] pnat, double[] v, double s, double bm1) {
        var ppr = new double[3] {0, 0, 0};
        S_Ab(pnat, v, s, bm1, ppr);
        return ppr;
    }
    [DllImport(DllFilename, EntryPoint = "iauAb", CallingConvention = CallingConvention.Cdecl)]
    static extern double S_Ab(double[] pnat, double[] v, double s, double bm1, [In, Out] double[] ppr);   
    
    
    /// <summary> Sofa.Apcg13(): For a GEOCENTRIC OBSERVER, prepare star−independent astrometry parameters
    /// for transformations between ICRS and GCRS coordinates. The caller supplies the date,
    /// and SOFA models are used to predict the Earth ephemeris. The parameters produced by this function
    /// are required in the parallax, light deflection and aberration parts of the astrometric
    /// transformation chain.</summary>
    /// <param name="date1">TDB as Julian date, part 1.
    /// TT maybe used without significant loss of accuracy.</param>
    /// <param name="date2">TDB as Julian date, part 2.
    /// TT maybe used without significant loss of accuracy.</param>
    /// <returns>Star−independent astrometry parameters, as a iauASTROM struct.
    /// Unchanged struct elements: along, xp1, yp1, sphi, cphi, diurab, era1, refa, refb.</returns>
    public static iauASTROM Apcg13(double date1, double date2) {
        var astrom = new iauASTROM();
        S_Apcg13(date1, date2, ref astrom);
        return astrom;
    }
    [DllImport(DllFilename, EntryPoint = "iauApcg13", CallingConvention = CallingConvention.Cdecl)]
    static extern void S_Apcg13(double date1, double date2, ref iauASTROM astrom); 
    
    
    /// <summary> Sofa.Apci13(): For a TERRESTRIAL OBSERVER, prepare star−independent astrometry parameters
    /// for transformations between ICRS and GCRS coordinates. The caller supplies the date,
    /// and SOFA models are used to predict the Earth ephemeris. The parameters produced by this function
    /// are required in the parallax, light deflection and aberration parts of the astrometric
    /// transformation chain.</summary>
    /// <param name="date1">TDB as Julian date, part 1.
    /// TT maybe used without significant loss of accuracy.</param>
    /// <param name="date2">TDB as Julian date, part 2.
    /// TT maybe used without significant loss of accuracy.</param>
    /// <returns> 2-Tuple:
    ///     astrom: Star−independent astrometry parameters, as a iauASTROM struct.
    ///             Unchanged struct elements: along, xp1, yp1, sphi, cphi, diurab, era1, refa, refb.
    ///     eo: equation of the origins (i.e., ERA-GST).</returns>
    public static Tuple<iauASTROM, double> Apci13(double date1, double date2) {
        var astrom = new iauASTROM();
        double eo = 0;
        S_Apci13(date1, date2, ref astrom, ref eo);
        var result = Tuple.Create(astrom, eo);
        return result;
    }
    [DllImport(DllFilename, EntryPoint = "iauApci13", CallingConvention = CallingConvention.Cdecl)]
    static extern void S_Apci13(double date1, double date2, ref iauASTROM astrom, ref double eo);


    /// <summary> Sofa.Apco13(): For a TERRESTRIAL OBSERVER, prepare star−independent astrometry parameters
    /// for transformations between ICRS and GCRS coordinates. The caller supplies UTC, site coordinates,
    /// ambient air conditions and observing wavelength, and SOFA models are used to predict
    /// the Earth ephemeris. The parameters produced by this function are required in the parallax,
    /// light deflection and aberration parts of the astrometric transformation chain.</summary>
    /// <param name="utc1">UTC as quasi-Julian date, part 1. Use Sofa.Dtf2d() to convert from calendar
    /// date and time of day into 2-part quasi-Julian date, as it implements the proper leap-second
    /// ambiguity convention.</param>
    /// <param name="utc2">UTC as quasi-Julian date, part 2.</param>
    /// <param name="dut1">UT1-UTC (seconds). Tabulated in IERS bulletins. </param>
    /// <param name="elong">Longitude of observer (radians, WGS84). CAUTION: east positive).</param>
    /// <param name="phi">Latitude of observer (radians, geodetic).</param>
    /// <param name="hm">Height above ellipsoid (meters, geodetic).</param>
    /// <param name="xp">Polar motion x coordinate (radians).
    /// Obtained from IERS bulletins, or may be set to zero with little consequence.</param>
    /// <param name="yp">Polar motion y coordinate (radians).
    /// Obtained from IERS bulletins, or may be set to zero with little consequence.</param>
    /// <param name="phpa">Pressure at the observer (hPa = mb) May be estimated from hm as
    /// 1013.25 * exp(-hm / (29.3 * tsl)) where tsl is approx. sea-level temperature Kelvin.</param>
    /// <param name="tc">Ambient temperature at the observer (deg C).</param>
    /// <param name="rh">Relative humidity at the observer (range 0-1).</param>
    /// <param name="w1">Wavelength of observation (micrometers).</param>
    /// <returns> 2-Tuple:
    ///     astrom: Star−independent astrometry parameters, as a iauASTROM struct.
    ///             Unchanged struct elements: along, xp1, yp1, sphi, cphi, diurab, era1, refa, refb.
    ///     eo: equation of the origins (i.e., ERA-GST).</returns>
    public static iauASTROM Apco13(double utc1, double utc2, double dut1, double elong,
        double phi, double hm, double xp, double yp, double phpa, double tc, double rh, double w1) {
        var astrom = new iauASTROM();
        S_Apco13(utc1, utc2, dut1, elong, phi, hm, xp, yp, phpa, tc, rh, w1, ref astrom);
        return astrom;
    }
    [DllImport(DllFilename, EntryPoint = "iauApco13", CallingConvention = CallingConvention.Cdecl)]
    static extern void S_Apco13(double utc1, double utc2, double dut1, double elong, double phi, double hm,
        double xp, double yp, double phpa, double tc, double rh, double w1, ref iauASTROM astrom);


    /// <summary> Sofa.Apcs13(): For an observer at known geocentric position and velocity,
    /// prepare star−independent astrometry parameters
    /// for transformations between ICRS and GCRS coordinates. The caller supplies the date,
    /// and SOFA models are used to predict the Earth ephemeris. The parameters produced by this function
    /// are required in the parallax, light deflection and aberration parts of the astrometric
    /// transformation chain.
    ///
    /// No assumptions are made about proximity to the earth, and this function may be used for deep space
    /// applications as well as Earth orbit and terrestrial locations.</summary>
    /// <param name="date1">TDB as Julian date, part 1.
    /// TT maybe used without significant loss of accuracy.</param>
    /// <param name="date2">TDB as Julian date, part 2.
    /// TT maybe used without significant loss of accuracy.</param>
    /// <param name="pv">Observer's geocentric position/velocity vectors, with respect to BCRS
    /// (Barycentric) axes, and NOT geocentric axes. (meters, meters/second).</param>
    /// <returns>Star−independent astrometry parameters, as a iauASTROM struct.
    /// Unchanged struct elements: along, xp1, yp1, sphi, cphi, diurab, era1, refa, refb.</returns>
    public static iauASTROM Apcs13(double date1, double date2, double[,] pv) {
        var astrom = new iauASTROM();
        S_Apcs13(date1, date2, pv, ref astrom);
        return astrom;
    }
    [DllImport(DllFilename, EntryPoint = "iauApcs13", CallingConvention = CallingConvention.Cdecl)]
    static extern void S_Apcs13(double date1, double date2, double[,] pv, ref iauASTROM astrom); 
    
        
    /// <summary> Sofa.Aper13(): For a TERRESTRIAL OBSERVER, use longitude to update "local" Earth rotation
    /// angle only.</summary>
    /// <param name="ut11">UT1 (not UTC) as Julian date, part 1.</param>
    /// <param name="ut12">UT1 (not UTC) as Julian date, part 2.</param>
    /// <param name="astrom">Astrometry struct iauASTROM, but only the along element is used.</param>
    /// <returns>astrom: Star−independent astrometry parameter era1, as the only updated element of
    /// this iauASTROM struct.</returns>
    public static iauASTROM Aper13(double ut11, double ut12, ref iauASTROM astrom) {
        S_Aper13(ut11, ut12, ref astrom);
        return astrom;
    }
    [DllImport(DllFilename, EntryPoint = "iauAper13", CallingConvention = CallingConvention.Cdecl)]
    static extern void S_Aper13(double ut11, double ut12, ref iauASTROM astrom);

    
    /// <summary> Sofa.Apio13(): For a TERRESTRIAL OBSERVER, prepare star−independent astrometry parameters
    /// for transformations between CIRS and observed coordinates. The caller supplies UTC, site coordinates,
    /// ambient air conditions and observing wavelength, and SOFA models are used to predict
    /// the Earth ephemeris. The parameters produced by this function are required in the parallax,
    /// light deflection and aberration parts of the astrometric transformation chain.</summary>
    /// <param name="utc1">UTC as quasi-Julian date, part 1. Use Sofa.Dtf2d() to convert from calendar
    /// date and time of day into 2-part quasi-Julian date, as it implements the proper leap-second
    /// ambiguity convention.</param>
    /// <param name="utc2">UTC as quasi-Julian date, part 2.</param>
    /// <param name="dut1">UT1-UTC (seconds). Tabulated in IERS bulletins. </param>
    /// <param name="elong">Longitude of observer (radians, WGS84). CAUTION: east positive).</param>
    /// <param name="phi">Latitude of observer (radians, geodetic).</param>
    /// <param name="hm">Height above ellipsoid (meters, geodetic).</param>
    /// <param name="xp">Polar motion x coordinate (radians).
    /// Obtained from IERS bulletins, or may be set to zero with little consequence.</param>
    /// <param name="yp">Polar motion y coordinate (radians).
    /// Obtained from IERS bulletins, or may be set to zero with little consequence.</param>
    /// <param name="phpa">Pressure at the observer (hPa = mb) May be estimated from hm as
    /// 1013.25 * exp(-hm / (29.3 * tsl)) where tsl is approx. sea-level temperature Kelvin.</param>
    /// <param name="tc">Ambient temperature at the observer (deg C).</param>
    /// <param name="rh">Relative humidity at the observer (range 0-1).</param>
    /// <param name="w1">Wavelength of observation (micrometers).</param>
    /// <returns> astrom: Star−independent astrometry parameters, as a iauASTROM struct.
    /// Struct elements unchanged: pmt, eb, eh, em, v, bm1, bpn.</returns>
    public static iauASTROM Apio13(double utc1, double utc2, double dut1, double elong,
        double phi, double hm, double xp, double yp, double phpa, double tc, double rh, double w1) {
        var astrom = new iauASTROM();
        S_Apio13(utc1, utc2, dut1, elong, phi, hm, xp, yp, phpa, tc, rh, w1, ref astrom);
        return astrom;
    }
    [DllImport(DllFilename, EntryPoint = "iauApio13", CallingConvention = CallingConvention.Cdecl)]
    static extern void S_Apio13(double utc1, double utc2, double dut1, double elong, double phi, double hm,
        double xp, double yp, double phpa, double tc, double rh, double w1, ref iauASTROM astrom);
    
    
    /// <summary>Sofa.Atcc13(): Transform a star’s ICRS catalog entry (epoch J2000.0) into ICRS
    ///                         astrometric place (update J2000.0 RA,Dec to datetime other than 2000.0).
    ///                         [NB: this fails for ATLAS refcat2, whose position epoch is NOT 2000.0.)
    ///                         </summary>
    /// <param name="rc">ICRS right ascension at J2000.0 (radians)</param>
    /// <param name="dc">ICRS declination at J2000.0 (radians)</param>
    /// <param name="pr">RA proper motion (radians/year; as dRA/dt, not cos(Dec)*dRA/dt)</param>
    /// <param name="pd">Dec proper motion (radians/year)</param>
    /// <param name="px">parallax (arcseconds)</param>
    /// <param name="rv">radial velocity (km/s, +ve if receding)</param>
    /// <param name="date1">Part 1 of 2-part Julian date.</param>
    /// <param name="date2">Part 2 of 2-part Julian date.</param>
    /// <returns>2-tuple of doubles: RA, Dec (ICRS astrometric).</returns>
    public static Tuple<double, double> Atcc13(double rc, double dc, double pr, double pd, 
        double px, double rv, double date1, double date2) {
        double ra = 0, da = 0;
        S_Atcc13(rc, dc, pr, pd, px, rv, date1, date2, ref ra, ref da);
        return Tuple.Create(ra, da); }
    [DllImport(DllFilename, EntryPoint = "iauAtcc13", CallingConvention = CallingConvention.Cdecl)]
    static extern double S_Atcc13(double rc, double dc, double pr, double pd, double px, double rv,
        double date1, double date2, ref double ra, ref double da);   

    
    /// <summary> Sofa.Atccq(): Quick transformation of a star’s ICRS catalog entry (epoch J2000.0)
    /// into ICRS astrometric place, given precomputed star−independent astrometry parameters.
    ///
    /// Use of this function is appropriate when efficiency is important and where many star positions
    /// are to be transformed for one date. The star−independent parameters can be obtained by calling
    /// one of the functions iauApci[13], iauApcg[13], iauApco[13] or iauApcs[13].</summary>
    /// <param name="rc">ICRS RA at J2000.0 (radians).</param>
    /// <param name="dc">ICRS Dec at J2000.0 (radians).</param>
    /// <param name="pr">RA proper motion (radians/year, in raw dRA/dt, not in cos(Dec)*dRA/dt).</param>
    /// <param name="pd">Dec proper motion (radians/year).</param>
    /// <param name="px">Parallax (arcseconds).</param>
    /// <param name="rv">Radial valocity (km/s, positive for receding).</param>
    /// <param name="astrom">Star-independent astrometry parameters, as a iauASTROM struct generated
    /// by another SOFA function..</param>
    /// <returns>2-tuple:  ra: ICRS astrometric Right Ascension (radians).
    ///                    dec: ICRS astrometric Declination (radians).</returns>
    public static Tuple<double, double> Atccq(double rc, double dc, double pr, double pd, 
        double px, double rv, iauASTROM astrom) {
        double ra = 0, da = 0;
        S_Atccq(rc, dc, pr, pd, px, rv, astrom, ref ra, ref da);
        return Tuple.Create(ra, da); }
    [DllImport(DllFilename, EntryPoint = "iauAtccq", CallingConvention = CallingConvention.Cdecl)]
    static extern double S_Atccq(double rc, double dc, double pr, double pd, double px, double rv,
        iauASTROM astrom, ref double ra, ref double da);   
    
    
    /// <summary>Sofa.Atci13(): Transform ICRS star data, epoch J2000.0, to CIRS.
    /// NB: SOFA param namings are reversed: rc & dc are ICRS, not CIRS; ri & di are CIRS.</summary>
    /// <param name="rc">ICRS right ascension at J2000.0 (radians)</param>
    /// <param name="dc">ICRS declination at J2000.0 (radians)</param>
    /// <param name="pr">RA proper motion (radians/year; as dRA/dt, not cos(Dec)*dRA/dt)</param>
    /// <param name="pd">Dec proper motion (radians/year)</param>
    /// <param name="px">parallax (arcseconds)</param>
    /// <param name="rv">radial velocity (km/s, +ve if receding)</param>
    /// <param name="date1">Part 1 of 2-part Julian date.</param>
    /// <param name="date2">Part 2 of 2-part Julian date.</param>
    /// <returns>3-tuple of doubles: RA, Dec (ICRS astrometric),
    ///                              eo (equation of origins, ERA-GST).</returns>
    public static Tuple<double, double, double> Atci13(double rc, double dc, double pr, double pd, 
        double px, double rv, double date1, double date2) {
        double ri = 0, di = 0, eo = 0;
        S_Atci13(rc, dc, pr, pd, px, rv, date1, date2, ref ri, ref di, ref eo);
        return Tuple.Create(ri, di, eo); }
    [DllImport(DllFilename, EntryPoint = "iauAtci13", CallingConvention = CallingConvention.Cdecl)]
    static extern double S_Atci13(double rc, double dc, double pr, double pd, double px, double rv,
        double date1, double date2, ref double ri, ref double di, ref double eo); 
    
    
    
    
    
    /// <summary>Sofa.Atco13(): ICRS RA,Dec to observed place. The caller supplies UTC,
    /// site coordinates, ambient air conditions and observing wavelength.</summary>
    /// <param name="rc">ICRS right ascension at J2000.0 (radians)</param>
    /// <param name="dc">ICRS declination at J2000.0 (radians)</param>
    /// <param name="pr">RA proper motion (radians/year; as dRA/dt, not cos(Dec)*dRA/dt).</param>
    /// <param name="pd">Dec proper motion (radians/year).</param>
    /// <param name="px">parallax (arcseconds).</param>
    /// <param name="rv">radial velocity (km/s, +ve if receding)</param>
    /// <param name="utc1">Part 1 of 2-part UTC quasi-Julian date.</param>
    /// <param name="utc2">Part 2 of 2-part UTC quasi-Julian date.</param>
    /// <param name="dut1">UT1−UTC (seconds, Note 5).</param>
    /// <param name="elong">Longitude (radians, WGS84, east +ve).</param>
    /// <param name="phi">Latitude (radians, geodetic, WGS84).</param>
    /// <param name="hm">Height above ellipsoid (meters, geodetic, WGS84).</param>
    /// <param name="xp">Polar motion coordinate (radians, from IERS bulletins, often ~zero).</param>
    /// <param name="yp">Polar motion coordinate (radians, from IERS bulletins, often ~zero).</param>
    /// <param name="phpa">Pressure at the observer (hPa = mB).</param>
    /// <param name="tc">Ambient temperature at the observer (degrees C).</param>
    /// <param name="rh">Relative humidity at the observer (in range [0-1]).</param>
    /// <param name="wl">Observation wavelength (micrometers).</param> 
    /// <returns>6-tuple of doubles: aob: observed azimuth (radians),
    ///                              zob: observed zenith distance NB: NOT altitude (radians),
    ///                              hob: observed hour angle (radians),
    ///                              dob: observed declination (radians),
    ///                              rob: observed right ascension (radians, CIO-based),
    ///                              eo:  equation of the origins (ERA, GST).</returns>
    public static Tuple<double, double, double, double, double, double> Atco13(double rc, double dc, 
        double pr, double pd, double px, double rv, double utc1, double utc2, double dut1, 
        double elong, double phi, double hm, double xp, double yp, double phpa, 
        double tc, double rh, double wl) {
        double aob = 0, zob = 0, hob = 0, dob = 0, rob = 0, eo = 0;
        S_Atco13(rc, dc, pr, pd, px, rv, utc1, utc2, dut1, elong, phi, hm, xp, yp, phpa, tc, rh, wl,
            ref aob, ref zob, ref hob, ref dob, ref rob, ref eo);
        return Tuple.Create(aob, zob, hob, dob, rob, eo); }
    [DllImport(DllFilename, EntryPoint = "iauAtco13", CallingConvention = CallingConvention.Cdecl)]
    static extern double S_Atco13(double rc, double dc, double pr, double pd, double px, double rv,
        double utc1, double utc2, double dut1, double elong, double phi, double hm, double xp, double yp,
        double phpa, double tc, double rh, double wl,
        ref double aob, ref double zob, ref double hob, ref double dob, ref double rob, ref double eo);
    
    
    /// <summary>Sofa.Atic13(): Transform star RA,Dec from geocentric CIRS to ICRS astrometric.</summary>
    /// <param name="ri">CIRS geocentric RA (radians)</param>
    /// <param name="di">CIRS geocentric Dec (radians)</param>
    /// <param name="date1">Part 1 of 2-part Julian date.</param>
    /// <param name="date2">Part 2 of 2-part Julian date.</param>
    /// <returns>3-tuple of doubles: rc, dc: RA, Dec (ICRS astrometric),
    ///                              eo (equation of origins, ERA-GST).</returns>
    public static Tuple<double, double, double> Atic13(double ri, double di, double date1, double date2) {
        double rc = 0, dc = 0, eo = 0;
        S_Atic13(ri, di, date1, date2, ref rc, ref dc, ref eo);
        return Tuple.Create(rc, dc, eo); }
    [DllImport(DllFilename, EntryPoint = "iauAtic13", CallingConvention = CallingConvention.Cdecl)]
    static extern double S_Atic13(double ri, double di, double date1, double date2, 
        ref double rc, ref double dc, ref double eo); 
    
    
    /// <summary>Sofa.Atio13(): Transform CIRS RA,Dec to observed place. The caller supplies UTC,
    ///                         site coordinates, ambient air conditions and observing wavelength.</summary>
    /// <param name="ri">CIRS geocentric RA (radians)</param>
    /// <param name="di">CIRS geocentric Dec (radians)</param>
    /// <param name="utc1">Part 1 of 2-part UTC quasi-Julian date.</param>
    /// <param name="utc2">Part 2 of 2-part UTC quasi-Julian date.</param>
    /// <param name="dut1">UT1−UTC (seconds, Note 5).</param>
    /// <param name="elong">Longitude (radians, WGS84, east +ve).</param>
    /// <param name="phi">Latitude (radians, geodetic, WGS84).</param>
    /// <param name="hm">Height above ellipsoid (meters, geodetic, WGS84).</param>
    /// <param name="xp">Polar motion coordinate (radians, from IERS bulletins, often ~zero).</param>
    /// <param name="yp">Polar motion coordinate (radians, from IERS bulletins, often ~zero).</param>
    /// <param name="phpa">Pressure at the observer (hPa = mB).</param>
    /// <param name="tc">Ambient temperature at the observer (degrees C).</param>
    /// <param name="rh">Relative humidity at the observer (in range [0-1]).</param>
    /// <param name="wl">Observation wavelength (micrometers).</param> 
    /// <returns>5-tuple of doubles: aob: observed azimuth (radians),
    ///                              zob: observed zenith distance (radians),
    ///                              hob: observed hour angle (radians),
    ///                              dob: observed declination (radians),
    ///                              rob: observed right ascension (radians, CIO-based)</returns>.
    public static Tuple<double, double, double, double, double> Atio13(double ri, double di, 
        double utc1, double utc2, double dut1, double elong, double phi, double hm, double xp, double yp,
        double phpa, double tc, double rh, double wl) {
        double aob = 0, zob = 0, hob = 0, dob = 0, rob = 0;
        S_Atio13(ri, di, utc1, utc2, dut1, elong, phi, hm, xp, yp, phpa, tc, rh, wl,
            ref aob, ref zob, ref hob, ref dob, ref rob);
        return Tuple.Create(aob, zob, hob, dob, rob); }
    [DllImport(DllFilename, EntryPoint = "iauAtio13", CallingConvention = CallingConvention.Cdecl)]
    static extern double S_Atio13(double ri, double di, double utc1, double utc2, 
        double dut1, double elong, double phi, double hm, double xp, double yp,
        double phpa, double tc, double rh, double wl,
        ref double aob, ref double zob, ref double hob, ref double dob, ref double rob); 
    
    
    /// <summary>Sofa.Atoc13(): Transform observed place at a ground-based site to ICRS
    ///                         astrometric RA,Dec. The caller supplies UTC, site coordinates,
    ///                         ambient air conditions, and observing wavelength.</summary>
    /// <param name="type">type of coordinates (string of length one, case-insensitive):
    ///                    "R" -> ob1=RA, ob2=Dec,
    ///                    "H" -> ob1=Hour Angle, ob2=Dec,
    ///                 or "A" -> ob1=Azimuth (N=0, E=90deg), ob2=zenith distance).</param>
    /// <param name="ob1">observed RA, Hour Angle, or Azimuth (N=0, E+) (radians)</param>
    /// <param name="ob2">observed Zenith Distance or Declination (radians)</param>
    /// <param name="utc1">Part 1 of 2-part UTC quasi-Julian date.</param>
    /// <param name="utc2">Part 2 of 2-part UTC quasi-Julian date.</param>
    /// <param name="dut1">UT1−UTC (seconds, Note 5).</param>
    /// <param name="elong">Longitude (radians, WGS84, east +ve).</param>
    /// <param name="phi">Latitude (radians, geodetic, WGS84).</param>
    /// <param name="hm">Height above ellipsoid (meters, geodetic, WGS84).</param>
    /// <param name="xp">Polar motion coordinate (radians, from IERS bulletins, often ~zero).</param>
    /// <param name="yp">Polar motion coordinate (radians, from IERS bulletins, often ~zero).</param>
    /// <param name="phpa">Pressure at the observer (hPa = mB).</param>
    /// <param name="tc">Ambient temperature at the observer (degrees C).</param>
    /// <param name="rh">Relative humidity at the observer (in range [0-1]).</param>
    /// <param name="wl">Observation wavelength (micrometers).</param> 
    /// <returns>2-tuple of doubles: rc, dc: RA, Dec (ICRS astrometric).</returns>.
    public static Tuple<double, double> Atoc13(string type, double ob1, double ob2, double utc1, double utc2, 
        double dut1, double elong, double phi, double hm, double xp, double yp, 
        double phpa, double tc, double rh, double wl) {
        double rc = 0, dc = 0;
        S_Atoc13(type, ob1, ob2, utc1, utc2, dut1, elong, phi, hm, xp, yp, phpa, tc, rh, wl,
            ref rc, ref dc);
        return Tuple.Create(rc, dc); }
    [DllImport(DllFilename, EntryPoint = "iauAtoc13", CallingConvention = CallingConvention.Cdecl)]
    static extern double S_Atoc13(string type, double ob1, double ob2, double utc1, double utc2, 
        double dut1, double elong, double phi, double hm, double xp, double yp, double phpa, double tc, 
        double rh, double wl, ref double rc, ref double dc); 
    
    
    /// <summary>Sofa.Atoi13(): Transform observed place to CIRS. The caller supplies UTC,
    ///                         site coordinates, ambient air conditions and observing wavelength.</summary>
    /// <param name="type">type of coordinates (string of length one, case-insensitive):
    ///                    "R" -> ob1=RA, ob2=Dec,
    ///                    "H" -> ob1=Hour Angle, ob2=Dec,
    ///                 or "A" -> ob1=Azimuth (N=0, E=90deg), ob2=zenith distance).</param>
    /// <param name="ob1">observed RA, Hour Angle, or Azimuth (N=0, E+) (radians)</param>
    /// <param name="ob2">observed Zenith Distance or Declination (radians)</param>
    /// <param name="utc1">Part 1 of 2-part UTC quasi-Julian date.</param>
    /// <param name="utc2">Part 2 of 2-part UTC quasi-Julian date.</param>
    /// <param name="dut1">UT1−UTC (seconds, Note 5).</param>
    /// <param name="elong">Longitude (radians, WGS84, east +ve).</param>
    /// <param name="phi">Latitude (radians, geodetic, WGS84).</param>
    /// <param name="hm">Height above ellipsoid (meters, geodetic, WGS84).</param>
    /// <param name="xp">Polar motion coordinate (radians, from IERS bulletins, often ~zero).</param>
    /// <param name="yp">Polar motion coordinate (radians, from IERS bulletins, often ~zero).</param>
    /// <param name="phpa">Pressure at the observer (hPa = mB).</param>
    /// <param name="tc">Ambient temperature at the observer (degrees C).</param>
    /// <param name="rh">Relative humidity at the observer (in range [0-1]).</param>
    /// <param name="wl">Observation wavelength (micrometers).</param> 
    /// <returns>2-tuple of doubles: ri, di: RA, Dec (radians, CIO-based).</returns>.
    public static Tuple<double, double> Atoi13(string type, double ob1, double ob2, double utc1, double utc2, 
        double dut1, double elong, double phi, double hm, double xp, double yp, 
        double phpa, double tc, double rh, double wl) {
        double ri = 0, di = 0;
        S_Atoi13(ref type, ob1, ob2, utc1, utc2, dut1, elong, phi, hm, xp, yp, phpa, tc, rh, wl,
            ref ri, ref di);
        return Tuple.Create(ri, di); }
    [DllImport(DllFilename, EntryPoint = "iauAtoi13", CallingConvention = CallingConvention.Cdecl)]
    static extern double S_Atoi13(ref string type, double ob1, double ob2, double utc1, double utc2, 
        double dut1, double elong, double phi, double hm, double xp, double yp, double phpa, double tc, 
        double rh, double wl, ref double ri, ref double di);
    
    
    /// <summary>Sofa.Pmsafe(): Star proper motion: update star catalog data for space motion, with
    ///                         special handling to handle the zero parallax case.</summary>
    /// <param name="ra1">Right Ascension, before (radians)</param>
    /// <param name="dec1">Declination, before (radians)</param>
    /// <param name="pmr1">RA Proper motion, before (radians/year).</param>
    /// <param name="pmd1">Declination Proper motion, before (radians/year).</param>
    /// <param name="px1">Parallax, before (arcseconds).</param>
    /// <param name="rv1">Radial velocity, before (km/s, +v -> receding).</param>
    /// <param name="ep1a">Part A of 2-part Epoch TDB quasi-JD, before.</param>
    /// <param name="ep1b">Part B of 2-part Epoch TDB quasi-JD, before.</param>
    /// <param name="ep2a">Part A of 2-part Epoch TDB quasi-JD, after.</param>
    /// <param name="ep2b">Part B of 2-part Epoch TDB quasi-JD, after.</param>
    /// <returns>6-tuple of doubles:  ra2:   Right Ascension, after (radians)
    ///                               dec2:  Declination, after (radians)
    ///                               pmr2:  RA proper motion, after (radians/year)
    ///                               pmd2:  Declination proper motion, after (radians/year)
    ///                               px2:   Parallax, after (arcseconds)
    ///                               rv2:   Radial velocity, after (km/s, +v -> receding). </returns>.
    public static Tuple<double, double, double, double, double, double> Pmsafe(
        double ra1, double dec1, double pmr1, double pmd1, double px1, double rv1, 
        double ep1a, double ep1b, double ep2a, double ep2b) {
        double ra2 = 0, dec2 = 0, pmr2 = 0, pmd2 = 0, px2 = 0, rv2 = 0;
        S_Pmsafe(ra1, dec1, pmr1, pmd1, px1, rv1, ep1a, ep1b, ep2a, ep2b, 
            ref ra2, ref dec2, ref pmr2, ref pmd2, ref px2, ref rv2);
        return Tuple.Create(ra2, dec2, pmr2, pmd2, px2, rv2); }
    [DllImport(DllFilename, EntryPoint = "iauPmsafe", CallingConvention = CallingConvention.Cdecl)]
    static extern double S_Pmsafe(double ra1, double dec1, double pmr1, double pmd1, 
        double px1, double rv1, double ep1a, double ep1b, double ep2a, double ep2b,
        ref double ra2, ref double dec2, ref double pmr2, ref double pmd2, ref double px2, ref double rv2);
    
    
    /// <summary>Sofa.Pvtob(): Position and velocity (CIRS) of a terrestrial observing station.</summary>
    /// <param name="elong">Longitude (radians, positive=east)</param>
    /// <param name="phi">Latitude (radians, geodetic)</param>
    /// <param name="hm">Height above reference ellipsoid (meters, geodetic).</param>
    /// <param name="xp">X-coordinate of the pole (radians, close to zero).</param>
    /// <param name="yp">Y-coordinate of the pole (radians, close to zero).</param>
    /// <param name="sp">TIO locator (radians, close to zero).</param>
    /// <param name="theta">Earth rotation angle (radians)</param>
    /// <returns>2-dim. array of doubles: position/velocity vector (m, m/s, CIRS).</returns>.
    public static double[,] Pvtob(double elong, double phi, double hm, double xp, double yp, 
        double sp, double theta) {
        var pv = new double[2,3] { {0, 0, 0}, {0, 0, 0} };
        S_Pvtob(elong, phi, hm, xp, yp, sp, theta, pv);
        return pv; }
    [DllImport(DllFilename, EntryPoint = "iauPvtob", CallingConvention = CallingConvention.Cdecl)]
    static extern double S_Pvtob(double elong, double phi, double hm, double xp, double yp, 
        double sp, double theta, [In, Out] double[,] pv);
    
    
    /// <summary>Sofa.Refco(): Determine the constants A and B in the atmospheric refraction model
    ///                        dZ = A tan Z + B tan^3 Z.</summary>
    /// <param name="phpa">Pressure at the observer (hPa = millibar)</param>
    /// <param name="tc">Ambient temperature at the observer (deg C)</param>
    /// <param name="rh">Relative humidity at the observer (range [0-1]).</param>
    /// <param name="wl">Wavelength of observation (micrometers).</param>
    /// <returns>2-tuple of doubles: refa: tan Z coefficient (radians)
    ///                              refb: tan^3 Z coefficient (radians).</returns>.
    public static Tuple<double, double> Refco(double phpa, double tc, double rh, double wl) {
        double refa = 0, refb = 0;
        S_Refco(phpa, tc, rh, wl, ref refa, ref refb);
        return Tuple.Create(refa, refb); }
    [DllImport(DllFilename, EntryPoint = "iauRefco", CallingConvention = CallingConvention.Cdecl)]
    static extern double S_Refco(double phpa, double tc, double rh, double wl, 
        ref double refa, ref double refb);

    #endregion    

#region "Astronomy/Ephemerides"   
    /*######### Astronomy/Ephemerides ####################################################################*/
    // No need to expose as public *for now* (reproduce these as needed, in C# within AstroLib):
    // .....................................................................................................

    
    /// <summary>Sofa.Epv00(): Earth position and velocity, heliocentric and barycentric, with
    ///                        respect to the Barycentric Celestial Reference System.</summary>
    /// <param name="date1">Part 1 of quasi-Julian form of TDB (or possibly TT) date.</param>
    /// <param name="date2">Part 2 of quasi-Julian form of TDB (or possibly TT) date.</param>
    /// <returns>2-tuple of 2-dim. arrays of doubles:
    ///              pvh[2][3]: heliocentric Earth position and velocity (AU, AU/day)
    ///              pvb[2][3]: barycentric Earth position and velocity (AU, AU/day) </returns>.
    public static Tuple<double[,], double[,]> Epv00(double date1, double date2) {
        var pvh = new double[2, 3] { {0, 0, 0}, {0, 0, 0} };
        var pvb = new double[2, 3] { {0, 0, 0}, {0, 0, 0} };
        var status = S_Epv00(date1, date2, pvh, pvb);
        return Tuple.Create(pvh, pvb); }
    [DllImport(DllFilename, EntryPoint = "iauEpv00", CallingConvention = CallingConvention.Cdecl)]
    static extern int S_Epv00(double date1, double date2, [In, Out] double[,] pvh, [In, Out] double[,] pvb); 
    
    
    /// <summary>Sofa.Moon98(): Approximate geocentric position and velocity of the Moon.</summary>
    /// <param name="date1">Part 1 of quasi-Julian form of TT (or possibly TDB) date.</param>
    /// <param name="date2">Part 2 of quasi-Julian form of TT (or possibly TDB) date.</param>
    /// <returns>2-dim. array of doubles:
    ///              pv[2][3]: Moon position and velocity, GCRS (in AU, AU/d)</returns>.
    public static double[,] Moon98(double date1, double date2) {
        var pv = new double[2, 3] { {0, 0, 0}, {0, 0, 0} };
        S_Moon98(date1, date2, pv);
        return pv; }
    [DllImport(DllFilename, EntryPoint = "iauMoon98", CallingConvention = CallingConvention.Cdecl)]
    static extern void S_Moon98(double date1, double date2, [In, Out] double[,] pv); 

    
    /// <summary>Sofa.Plan94(): Approximate heliocentric position and velocity of a nominated major
    ///                         planet: Mercury, Venus, EMB, Mars, Jupiter, Saturn, Uranus or
    ///                         Neptune (but not the Earth itself).</summary>
    /// <param name="date1">Part 1 of quasi-Julian form of TDB (or possibly TT) date.</param>
    /// <param name="date2">Part 2 of quasi-Julian form of TDB (or possibly TT) date.</param>
    /// <param name="np">planet (1=Mercury, 2=Venus, 3=EMB (Earth-Moon Barycenter; for earth alone use
    ///                         Epv00), 4=Mars, 5=Jupiter, 6=Saturn, 7=Uranus, 8=Neptune)</param>
    /// <returns>2-dim. array of doubles:
    ///              pv[2][3]: Planet position and velocity, heliocentric J2000.0 (in AU, AU/d)</returns>.
    public static double[,] Plan94(double date1, double date2, int np) {
        var pv = new double[2, 3] { {0, 0, 0}, {0, 0, 0} };
        var status = S_Plan94(date1, date2, np, pv);
        return pv; }
    [DllImport(DllFilename, EntryPoint = "iauPlan94", CallingConvention = CallingConvention.Cdecl)]
    static extern int S_Plan94(double date1, double date2, int np, [In, Out] double[,] pv);
#endregion    

#region "Astronomy/FundamentalArgs"
    /*######### Astronomy/FundamentalArgs ################################################################*/
    // No need to expose as public *for now* (reproduce these as needed, in C# within AstroLib):
    //      iauFad03  Fundamental argument, IERS Conventions (2003): mean elongation of the Moon
    //                    from the Sun.
    //      iauFae03  Fundamental argument, IERS Conventions (2003): mean longitude of Earth.
    //      iauFaf03  Fundamental argument, IERS Conventions (2003): mean longitude of the Moon
    //                    minus mean longitude of the ascending node.
    //      iauFaju03 Fundamental argument, IERS Conventions (2003): mean longitude of Jupiter.
    //      iauFal03  Fundamental argument, IERS Conventions (2003): mean anomaly of the Moon.
    //      iauFalp03 Fundamental argument, IERS Conventions (2003): mean anomaly of the Sun.
    //      iauFama03 Fundamental argument, IERS Conventions (2003): mean longitude of Mars.
    //      iauFame03 Fundamental argument, IERS Conventions (2003): mean longitude of Mercury.
    //      iauFane03 Fundamental argument, IERS Conventions (2003): mean longitude of Neptune.
    //      iauFaom03 Fundamental argument, IERS Conventions (2003): mean longitude of the
    //                    Moon’s ascending node.
    //      iauFapa03 Fundamental argument, IERS Conventions (2003): general accumulated
    //                    precession in longitude.
    //      iauFasa03 Fundamental argument, IERS Conventions (2003): mean longitude of Saturn.
    //      iauFaur03 Fundamental argument, IERS Conventions (2003): mean longitude of Uranus.
    //      iauFave03 Fundamental argument, IERS Conventions (2003): mean longitude of Venus.
    //......................................................................................................
#endregion    
    
#region "Astronomy/PrecNutPolar"
    /*######### Astronomy/PrecNutPolar ###################################################################*/
    // No need to expose as public *for now* (reproduce these as needed, in C# within AstroLib):
    //      iauAtoiq  Quick observed place to CIRS, given the star−independent astrometry parameters.
    //      iauBi00   Frame bias components of IAU 2000 precession−nutation models; part of the
    //                    Mathews−Herring−Buffett (MHB2000) nutation series, with additions.
    //      iauBp00   Frame bias and precession, IAU 2000.
    //      iauBp06   Frame bias and precession, IAU 2006.
    //      iauBpn2xy Extract from the bias−precession−nutation matrix the X,Y coordinates of
    //                    the Celestial Intermediate Pole.
    //      iauC2i00a Form the celestial−to−intermediate matrix for a given date using the
    //                    IAU 2000A precession−nutation model.
    //      iauC2i00b Form the celestial−to−intermediate matrix for a given date using the
    //                    IAU 2000B precession−nutation model.
    //      iauC2i06a Form the celestial−to−intermediate matrix for a given date using the
    //                    IAU 2006 precession and IAU 2000A nutation models.
    //      iauC2ibpn Form the celestial−to−intermediate matrix for a given date given the
    //                    bias−precession−nutation matrix. IAU 2000.
    //      iauC2ixy  Form the celestial to intermediate−frame−of−date matrix for a given date
    //                    when the CIP X,Y coordinates are known. IAU 2000.
    //      iauC2ixys Form the celestial to intermediate−frame−of−date matrix given the CIP X,Y
    //                    and the CIO locator s.
    //      iauC2t00a Form the celestial to terrestrial matrix given the date, the UT1 and the
    //                    polar motion, using the IAU 2000A precession−nutation model.
    //      iauC2t00b Form the celestial to terrestrial matrix given the date, the UT1 and the
    //                    polar motion, using the IAU 2000B precession−nutation model.
    //      iauC2t06a Form the celestial to terrestrial matrix given the date, the UT1 and the
    //                    polar motion, using the IAU 2006/2000A precession−nutation model.
    //      iauC2tcio Assemble the celestial to terrestrial matrix from CIO−based components
    //                    (the celestial−to−intermediate matrix, the Earth Rotation Angle and
    //                    the polar motion matrix).
    //      iauC2teqx Assemble the celestial to terrestrial matrix from equinox−based components
    //                    (the celestial−to−true matrix, the Greenwich Apparent Sidereal Time and
    //                    the polar motion matrix).
    //      iauC2tpe  Form the celestial to terrestrial matrix given the date, the UT1, the nutation
    //                    and the polar motion. IAU 2000.
    //      iauC2txy  Form the celestial to terrestrial matrix given the date, the UT1,
    //                    the CIP coordinates and the polar motion. IAU 2000.
    //      iauEo06a  Equation of the origins, IAU 2006 precession and IAU 2000A nutation.
    //      iauEors   Equation of the origins, given the classical NPB matrix and the quantity s.
    //      iauFw2m   Form rotation matrix given the Fukushima−Williams angles.
    //      iauFw2xy  CIP X,Y given Fukushima−Williams bias−precession−nutation angles.
    //      iauLtp    Long−term (ca. 40K year) precession matrix.
    //      iauLtpb   Long−term (ca. 40K year) precession matrix, including ICRS frame bias.
    //      iauLtpecl Long−term (ca. 40K year) precession of the ecliptic.
    //      iauLtpequ Long−term (ca. 40K year) precession of the equator.
    //      iauNum00a Form the matrix of nutation for a given date, IAU 2000A model.
    //      iauNum00b Form the matrix of nutation for a given date, IAU 2000B model.
    //      iauNum06a Form the matrix of nutation for a given date, IAU 2006/2000A model.
    //      iauNumat  Form the matrix of nutation.
    //      iauNut00a Nutation, IAU 2000A model (MHB2000 luni−solar and planetary nutation
    //                    with free core nutation omitted).
    //      iauNut00b Nutation, IAU 2000B model.
    //      iauNut06a IAU 2000A nutation with adjustments to match the IAU 2006 precession.
    //      iauNut80  Nutation, IAU 1980 model.
    //      iauNutm80 Form the matrix of nutation for a given date, IAU 1980 model.
    //      iauObl06  Mean obliquity of the ecliptic, IAU 2006 precession model.
    //      iauObl80  Mean obliquity of the ecliptic, IAU 1980 precession model.
    //      iauP06e   Precession angles, IAU 2006, equinox based.
    //      iauPb06   This function forms three Euler angles which implement general precession
    //                    from epoch J2000.0, using the IAU 2006 model. Frame bias (the offset
    //                    between ICRS and mean J2000.0) is included.
    //      iauPfw06  Precession angles, IAU 2006 (Fukushima−Williams 4−angle formulation).
    //      iauPmat00 Precession matrix (including frame bias) from GCRS to a specified date,
    //                    IAU 2000 model.
    //      iauPmat06 Precession matrix (including frame bias) from GCRS to a specified date,
    //                    IAU 2006 model.
    //      iauPmat76 Precession matrix from J2000.0 to a specified date, IAU 1976 model.
    //      iauPn00   Precession−nutation, IAU 2000 model: a multi−purpose function, supporting
    //                    classical (equinox−based) use directly and CIO−based use indirectly.
    //      iauPn00a  Precession−nutation, IAU 2000A model: a multi−purpose function, supporting
    //                    classical (equinox−based) use directly and CIO−based use indirectly.
    //      iauPn00b  Precession−nutation, IAU 2000B model: a multi−purpose function, supporting
    //                    classical (equinox−based) use directly and CIO−based use indirectly.
    //      iauPn06   Precession−nutation, IAU 2006 model: a multi−purpose function, supporting
    //                    classical (equinox−based) use directly and CIO−based use indirectly.
    //      iauPn06a  Precession−nutation, IAU 2006/2000A models: a multi−purpose function, supporting
    //                    classical (equinox−based) use directly and CIO−based use indirectly.
    //      iauPnm00a Form the matrix of precession−nutation for a given date (including frame bias),
    //                    equinox based, IAU 2000A model.
    //      iauPnm00b Form the matrix of precession−nutation for a given date (including frame bias),
    //                    equinox based, IAU 2000B model.
    //      iauPnm06a Form the matrix of precession−nutation for a given date (including frame bias),
    //                    equinox based, IAU 2006 precession and IAU 2000A nutation models.
    //      iauPnm80  Form the matrix of precession/nutation for a given date, IAU 1976 precession model,
    //                    IAU 1980 nutation model.
    //      iauPom00  Form the matrix of polar motion for a given date, IAU 2000.
    //      iauPr00   Precession−rate part of the IAU 2000 precession−nutation models (part of MHB2000).
    //      iauPrec76 IAU 1976 precession model.
    //      iauS00    The CIO locator s, positioning the Celestial Intermediate Origin on the equator
    //                    of the Celestial Intermediate Pole, given the CIP’s X,Y coordinates.
    //                    Compatible with IAU 2000A precession−nutation.
    //      iauS00a   The CIO locator s, positioning the Celestial Intermediate Origin on the equator
    //                    of the Celestial Intermediate Pole, using the IAU 2000A precession−nutation model.
    //      iauS00b   The CIO locator s, positioning the Celestial Intermediate Origin on the equator
    //                    of the Celestial Intermediate Pole, using the IAU 2000B precession−nutation model.
    //      iauS06    The CIO locator s, positioning the Celestial Intermediate Origin on the equator
    //                    of the Celestial Intermediate Pole, given the CIP’s X,Y coordinates.
    //                    Compatible with IAU 2006/2000A precession−nutation.
    //      iauS06a   The CIO locator s, positioning the Celestial Intermediate Origin on the equator
    //                    of the Celestial Intermediate Pole, using the IAU 2006 precession and
    //                    IAU 2000A nutation models.
    //      iauSp00   The TIO locator s’, positioning the Terrestrial Intermediate Origin on the equator
    //                    of the Celestial Intermediate Pole.
    //      iauXy06   X,Y coordinates of celestial intermediate pole from series based on
    //                    IAU 2006 precession and IAU 2000A nutation.
    //      iauXys00a For a given TT date, compute the X,Y coordinates of the Celestial Intermediate Pole
    //                    and the CIO locator s, using the IAU 2000A precession−nutation model.
    //      iauXys00b For a given TT date, compute the X,Y coordinates of the Celestial Intermediate Pole
    //                    and the CIO locator s, using the IAU 2000B precession−nutation model.
    //      iauXys06a For a given TT date, compute the X,Y coordinates of the Celestial Intermediate Pole
    //                    and the CIO locator s, using the IAU 2006 precession and
    //                    IAU 2000A nutation models.
    // .....................................................................................................
#endregion
    
#region "Astronomy/RotationAndTime"
    /*######### Astronomy/RotationAndTime ################################################################*/
    // No need to expose as public *for now* (reproduce these as needed, in C# within AstroLib):
    //      iauEe00   The equation of the equinoxes, compatible with IAU 2000 resolutions, given
    //                    the nutation in longitude and the mean obliquity.
    //      iauEe00a  Equation of the equinoxes, compatible with IAU 2000 resolutions.
    //      iauEe00b  Equation of the equinoxes, compatible with IAU 2000 resolutions but using
    //                    the truncated nutation model IAU 2000B.
    //      iauEe06a  Equation of the equinoxes, compatible with IAU 2000 resolutions and
    //                    IAU 2006/2000A precession−nutation.
    //      iauEect00 Equation of the equinoxes complementary terms, consistent with IAU 2000 resolutions.
    //      iauEqeq94 Equation of the equinoxes, IAU 1994 model.
    //      iauGmst00 Greenwich mean sidereal time (model consistent with IAU 2000 resolutions).
    //      iauGmst82 Universal Time to Greenwich mean sidereal time (IAU 1982 model).
    //      iauGst00a Greenwich apparent sidereal time (consistent with IAU 2000 resolutions).
    //      iauGst00b Greenwich apparent sidereal time (consistent with IAU 2000resolutions but
    //                    using the truncated nutation model IAU 2000B).
    //      iauGst06  Greenwich apparent sidereal time, IAU 2006, given the NPB matrix.
    //      iauGst94  Greenwich apparent sidereal time (consistent with IAU 1982/94 resolutions).
    // .....................................................................................................
    
    
    /// <summary>Sofa.Era00(): Earth rotation angle (IAU 2000 model).</summary>
    /// <param name="date1">Part 1 of Julian Date, UTC.</param>
    /// <param name="date2">Part 2 of Julian Date, UTC.</param>
    /// <returns>Earth Rotation Angle (radians)</returns>.
    public static double Era00(double date1, double date2) {
        var era = S_Era00(date1, date2);
        return era; }
    [DllImport(DllFilename, EntryPoint = "iauEra00", CallingConvention = CallingConvention.Cdecl)]
    static extern double S_Era00(double date1, double date2); 

    
    /// <summary>Sofa.Gmst06(): Greenwich mean sidereal time
    ///                         (consistent with IAU 2006 precession).
    /// Per SOFA manual "...if UT1 is used for both purposes, errors of order 100 microarcseconds result",
    ///     and such errors are almost always acceptable, we preserve herein the option to provide both.
    /// </summary>
    /// <param name="uta">Part 1 of Julian Date, UT1.</param>
    /// <param name="utb">Part 2 of Julian Date, UT1.</param>
    /// <param name="tta">Part 1 of Julian Date, TT (or use UT1 for milliarcsecond accuracy).</param>
    /// <param name="ttb">Part 2 of Julian Date, TT (or use UT1 for milliarcsecond accuracy).</param>
    /// <returns>Greenwich mean sidereal time (radians)</returns>.
    public static double Gmst06(double uta, double utb, double tta, double ttb) {
        var gmst = S_Gmst06(uta, utb, tta, ttb);
        return gmst; }
    [DllImport(DllFilename, EntryPoint = "iauGmst06", CallingConvention = CallingConvention.Cdecl)]
    static extern double S_Gmst06(double uta, double utb, double tta, double ttb); 
    
    
    /// <summary>Sofa.Gst06a(): Greenwich apparent sidereal time
    ///                         (consistent with IAU 2000 and 2006 resolutions).</summary>
    /// <param name="uta">Part 1 of Julian Date, UT1.</param>
    /// <param name="utb">Part 2 of Julian Date, UT1.</param>
    /// <param name="tta">Part 1 of Julian Date, TT (or use UT1 for milliarcsecond accuracy).</param>
    /// <param name="ttb">Part 2 of Julian Date, TT (or use UT1 for milliarcsecond accuracy).</param>
    /// <returns>Greenwich mean sidereal time (radians)</returns>.
    public static double Gst06a(double uta, double utb, double tta, double ttb) {
        var gast = S_Gst06a(uta, utb, tta, ttb);
        return gast; }
    [DllImport(DllFilename, EntryPoint = "iauGst06a", CallingConvention = CallingConvention.Cdecl)]
    static extern double S_Gst06a(double uta, double utb, double tta, double ttb);

    #endregion    
    
#region "Astronomy/SpaceMotion"
    /*######### Astronomy/SpaceMotion ####################################################################*/
    // No need to expose as public *for now* (reproduce these as needed, in C# within AstroLib):
    //      iauPvstar Convert star position+velocity vector to catalog coordinates.
    //      iauStarpv Convert star catalog coordinates to position+velocity vector.
    // .....................................................................................................
    
#endregion    
    
#region "Astronomy/StarCatalogs"
    /*######### Astronomy/StarCatalogs ###################################################################*/
    // No need to expose as public *for now* (reproduce these as needed, in C# within AstroLib):
    //      iauFk425  Convert B1950.0 FK4 star catalog data to J2000.0 FK5.
    //      iauFk45z  Convert a B1950.0 FK4 star position to J2000.0 FK5, assuming zero proper motion
    //                    in the FK5 system.
    //      iauFk524  Convert J2000.0 FK5 star catalog data to B1950.0 FK4.
    //      iauFk52h  Transform FK5 (J2000.0) star data into the Hipparcos system.
    //      iauFk54z  Convert a J2000.0 FK5 star position to B1950.0 FK4, assuming zero proper motion
    //                    in FK5 and parallax.
    //      iauFk5hip FK5 to Hipparcos rotation and spin.
    //      iauFk5hz  Transform an FK5 (J2000.0) star position into the system of the Hipparcos catalogue,
    //                    assuming zero Hipparcos proper motion.
    //      iauH2fk5  Transform Hipparcos star data into the FK5 (J2000.0) system.
    //      iauHfk5z  Transform a Hipparcos star position into FK5 J2000.0, assuming
    //                    zero Hipparcos proper motion.
    // .....................................................................................................
    
    
    /// <summary>Sofa.Starpm(): Star proper motion: update (for new time) star catalog data
    ///                         for space motion.
    /// NB: not safe for zero parallax and/or radial velocity: for these, always use Sofa.Pmsafe().
    /// </summary>
    /// <param name="ra1">Right Ascension, before (radians).</param>
    /// <param name="dec1">Declination, before (radians).</param>
    /// <param name="pmr1">RA proper motion, before (radians/year).</param>
    /// <param name="pmd1">Declination proper motion, before (radians/year).</param>
    /// <param name="px1">Parallax, before (arcseconds).</param>
    /// <param name="rv1">Radial velocity, before (km/s, positive=receding).</param>
    /// <param name="ep1a">Part 1 of Julian Date, before (TDB).</param>
    /// <param name="ep1b">Part 2 of Julian Date, before (TDB).</param>
    /// <param name="ep2a">Part 1 of Julian Date, after (TDB).</param>
    /// <param name="ep2b">Part 2 of Julian Date, after (TDB).</param>
    /// <returns>6-tuple of doubles: ra2:  Right Ascension, after (radians)
    ///                              dec2: Declination, after (radians)
    ///                              pmr2: RA proper motion, after (radians/year)
    ///                              pmd2: Declination proper motion, after (radians/year)
    ///                              px2:  Parallax, after (arcseconds)
    ///                              rv2:  Radial velocity, after (km/s, positive=receding).</returns>.
    public static Tuple<double, double, double, double, double, double> Starpm(double ra1, double dec1, 
        double pmr1, double pmd1, double px1, double rv1,
        double ep1a, double ep1b, double ep2a, double ep2b) {
        if (px1 == 0.0 && rv1 == 0.0) {
            throw new ArgumentException(
                "Use Sofa.Pmsafe() for zero parallax and/or zero radial velocity.");
        }
        double ra2 = 0, dec2 = 0, pmr2 = 0, pmd2 = 0, px2 = 0, rv2 = 0;
        var status = S_Starpm(ra1, dec1, pmr1, pmd1, px1, rv1, ep1a, ep1b, ep2a, ep2b,
            ref ra2, ref dec2,  ref pmr2, ref pmd2, ref px2, ref rv2);
        return Tuple.Create(ra2, dec2, pmr2, pmd2, px2, rv2); }
    [DllImport(DllFilename, EntryPoint = "iauStarpm", CallingConvention = CallingConvention.Cdecl)]
    static extern int S_Starpm(double ra1, double dec1, double pmr1, double pmd1, double px1, double rv1,
        double ep1a, double ep1b, double ep2a, double ep2b, 
        ref double ra2, ref double dec2, ref double pmr2, ref double pmd2, ref double px2, ref double rv2);

    #endregion    
    
#region "Astronomy/EclipticCoordinates"
    /*######### Astronomy/EclipticCoordinates ############################################################*/
    // No need to expose as public *for now* (reproduce these as needed, in C# within AstroLib):
    //      
    //      iauEcm06  ICRS equatorial to ecliptic rotation matrix, IAU 2006.
    //      iauLteceq Transformation from ecliptic coordinates (mean equinox and ecliptic of date)
    //                    to ICRS RA,Dec, using a long−term (ca. 40K year) precession model.
    //      iauLtecm  ICRS equatorial to ecliptic rotation matrix, long−term (ca. 40K year).
    //      iauLteqec Transformation from ICRS equatorial coordinates to ecliptic coordinates
    //                    (mean equinox and ecliptic of date) using a long−term (ca. 40K year)
    //                    precession model.
    // .....................................................................................................
    
    
    /// <summary>Sofa.Eceq06(): Transformation from ecliptic coordinates (mean equinox and ecliptic
    ///                         of date) to ICRS RA,Dec, using the IAU 2006 precession model.</summary>
    /// <param name="date1">Part 1 of 2-Part TT quasi-Julian date</param>
    /// <param name="date2">Part 2 of 2-Part TT quasi-Julian date</param>
    /// <param name="dl">Ecliptic longitude (radians)</param>
    /// <param name="db">Ecliptic latitude (radians)</param>
    /// <returns>2-tuple of doubles: dr, ICRS Right Ascension (radians)
    ///                              dd, ICRS Declination (radians)</returns>
    public static Tuple<double, double> Eceq06(double date1, double date2, double dl, double db) {
        double dr = 0, dd = 0;
        var status = S_Eceq06(date1, date2, dl, db, ref dr, ref dd);
        return Tuple.Create(dr, dd); }
    [DllImport(DllFilename, EntryPoint = "iauEceq06", CallingConvention = CallingConvention.Cdecl)]
    static extern int S_Eceq06(double date1, double date2, double dl, double db,
        ref double dr, ref double dd);
    
    
    /// <summary>Sofa.Eqec06(): Transformation from ICRS equatorial coordinates to ecliptic coordinates
    ///                         (mean equinox and ecliptic of date), using the IAU 2006 precession model.
    /// </summary>
    /// <param name="date1">Part 1 of 2-Part TT quasi-Julian date</param>
    /// <param name="date2">Part 2 of 2-Part TT quasi-Julian date</param>
    /// <param name="dr">ICRS Right Ascension (radians)</param>
    /// <param name="dd">ICRS Declination (radians)</param>
    /// <returns>2-tuple of doubles: dl, Ecliptic Longitude (radians)
    ///                              db, Ecliptic latitude (radians)</returns>
    public static Tuple<double, double> Eqec06(double date1, double date2, double dr, double dd) {
        double dl = 0, db = 0;
        S_Eqec06(date1, date2, dr, dd, ref dl, ref db);
        return Tuple.Create(dl, db); }
    [DllImport(DllFilename, EntryPoint = "iauEqec06", CallingConvention = CallingConvention.Cdecl)]
    static extern void S_Eqec06(double date1, double date2, double dr, double dd,
        ref double dl, ref double db);
#endregion    

#region "Astronomy/GalacticCoordinates"
    /*######### Astronomy/GalacticCoordinates ############################################################*/
    // No need to expose as public *for now* (reproduce these as needed, in C# within AstroLib):
    //      
    // .....................................................................................................
    
    
    /// <summary>Sofa.Icrs2g(): Transformation from ICRS to Galactic Coordinates.</summary>
    /// <param name="dr">ICRS right ascension (radians)</param>
    /// <param name="dd">ICRS Declination (radians)</param>
    /// <returns>2-tuple of doubles:
    ///              dl: Galactic longitude (radians).
    ///              db: Galactic latitude (radians).</returns>.
    public static Tuple<double, double> Icrs2g(double dr, double dd) {
        double dl = 0, db = 0;
        S_Icrs2g(dr, dd, ref dl, ref db);
        return Tuple.Create(dl, db); }
    [DllImport(DllFilename, EntryPoint = "iauIcrs2g", CallingConvention = CallingConvention.Cdecl)]
    static extern void S_Icrs2g(double dr, double dd, ref double dl, ref double db); 
    
    
    /// <summary>Sofa.G2icrs(): Transformation from Galactic Coordinates to ICRS.</summary>
    /// <param name="dl">Galactic longitude (radians).</param>
    /// <param name="db">Galactic latitude (radians).</param>
    /// <returns>2-tuple of doubles:
    ///              dr: ICRS right ascension (radians)
    ///              dd: ICRS Declination (radians) </returns>.
    public static Tuple<double, double> G2icrs(double dl, double db) {
        double dr = 0, dd = 0;
        S_G2icrs(dl, db, ref dr, ref dd);
        return Tuple.Create(dr, dd); }
    [DllImport(DllFilename, EntryPoint = "iauG2icrs", CallingConvention = CallingConvention.Cdecl)]
    static extern void S_G2icrs(double dl, double db, ref double dr, ref double dd); 

    
    
#endregion    
    
#region "Astronomy/GeodeticGeocentric"
    /*######### Astronomy/GeodeticGeocentric #############################################################*/
    // No need to expose as public *for now* (reproduce these as needed, in C# within AstroLib):
    //      iauEform  Earth reference ellipsoids.
    //      iauGc2gde Transform geocentric coordinates to geodetic for a reference ellipsoid
    //                    of specified form.
    //      iauGd2gce Transform geodetic coordinates to geocentric for a reference ellipsoid
    //                    of specified form.
    // .....................................................................................................
    
    
    /// <summary>Sofa.Gc2gd(): Transform geocentric coordinates to geodetic using the specified
    ///                        reference ellipsoid.
    ///                        NB: almost always use n=1 for WGS84.</summary>
    /// <param name="n">Ellipsoid identifier; n=1 -> WGS84, the normal case.</param>
    /// <param name="xyz">Geocentric vector (all 3 values in meters).</param>
    /// <returns>3-tuple of doubles:
    ///              elong: Longitude (radians)
    ///              phi: Latitude (radians)
    ///              height: Height above ellipsoid (meters, geodetic). </returns>
    public static Tuple<double, double, double> Gc2gd(int n, double[] xyz) {
        double elong = 0, phi = 0, height = 0;
        var status = S_Gc2gd(n, xyz, ref elong, ref phi, ref height);
        return Tuple.Create(elong, phi, height); }
    [DllImport(DllFilename, EntryPoint = "iauGc2gd", CallingConvention = CallingConvention.Cdecl)]
    static extern int S_Gc2gd(int n, double[] xyz, ref double elong, ref double phi, ref double height); 
    
    
    /// <summary>Sofa.Gd2gc(): Transform geodetic coordinates to geocentric using the specified
    ///                        reference ellipsoid.</summary>
    /// <param name="n">Ellipsoid identifier; n=1 -> WGS84, the normal case.</param>
    /// <param name="elong">Longitude (radians).</param>
    /// <param name="phi">Latitude (radians).</param>
    /// <param name="height">Height above ellipsoid (meters, geodetic).</param>
    /// <returns>Geocentric vector (all 3 values in meters).</returns>
    public static double[] Gd2gc(int n, double elong, double phi, double height) {
        var xyz = new double[3] {0, 0, 0}; 
        var status = S_Gd2gc(n, elong, phi, height, xyz);
        return xyz; }
    [DllImport(DllFilename, EntryPoint = "iauGd2gc", CallingConvention = CallingConvention.Cdecl)]
    static extern int S_Gd2gc(int n, double elong, double phi, double height, [In, Out] double[] xyz); 

#endregion    
    
#region "Astronomy/Timescales"
    /*######### Astronomy/Timescales #####################################################################*/
    // No need to expose as public *for now* (reproduce these as needed, in C# within AstroLib):
    //      iauD2tf   Decompose days to hours, minutes, seconds, fraction.
    //      iauDtdb   Time scale transformation:  Barycentric Dynamical Time, TDB, to Terrestrial Time, TT.
    //                     (Unclear or erroneous documentation on param "dtr") (~2 millisecond difference).
    //      iauTdbtt  Time scale transformation: Barycentric Dynamical Time, TDB, to Terrestrial Time, TT.
    //                    (Unclear or erroneous documentation on param "dtr") (~2 millisecond difference).
    //      iauTttdb  Time scale transformation: Terrestrial Time, TT, to Barycentric Dynamical Time, TDB.
    //                    (Unclear or erroneous documentation on param "dtr") (~2 millisecond difference).
    // .....................................................................................................
    
    
    /// <summary>Sofa.D2dtf(): Format for output a 2-part Julian Date
    /// (or in the case of UTC a quasi-JD form that includes special provision for leap seconds).</summary>
    /// <param name="scale">Time scale ID "UTC"</param>
    /// <param name="ndp">resolution</param>
    /// <param name="d1">first part of time as a 2-part Julian Date</param>
    /// <param name="d2">second part of time as a 2-part Julian Date</param>
    /// <returns>4-tuple: (int year, int month, int day,
    ///                    int array [hours, minutes, seconds, fraction])</returns>
    public static Tuple<int, int, int, int[]> D2dtf(string scale, int ndp, double d1, double d2) {
        int iy = 0, im = 0, id = 0;
        var ihmsf = new Int32[4] {0, 0, 0, 0};
        var status = S_D2dtf(scale, ndp, d1, d2, ref iy, ref im, ref id, ihmsf);
        return Tuple.Create(iy, im, id, ihmsf); }
    [DllImport(DllFilename, EntryPoint = "iauD2dtf", CallingConvention = CallingConvention.Cdecl)]
    static extern int S_D2dtf(string scale, int ndp, double d1, double d2, 
        ref int iy, ref int im, ref int id, [In, Out] int[] ihmsf);
    
    
    /// <summary>Sofa.Dat(): For a given UTC date, calculate Delta(AT) = TAI−UTC</summary>
    /// <param name="iy">Year (UTC)"</param>
    /// <param name="im">Month (UTC)</param>
    /// <param name="id">Day (UTC)</param>
    /// <param name="fd">Fraction of day (UTC)</param>
    /// <returns>double: TAI minus UTC (seconds)</returns>
    //
    // From the SOFA release 18 Manual, iauDat entry:
    // :−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−:
    // :
    // :                IMPORTANT
    // :
    // : A new version of this function must be
    // : produced whenever a new leap second is
    // : announced. There are four items to
    // : change on each such occasion:
    // :
    // : 1) A new line must be added to the set
    // : of statements that initialize the
    // : array "changes".
    // :
    // : 2) The constant IYV must be set to the
    // : current year.
    // :
    // : 3) The "Latest leap second" comment
    // : below must be set to the new leap
    // : second date.
    // :
    // : 4) The "This revision" comment, later,
    // : must be set to the current date.
    // :
    // : Change (2) must also be carried out
    // : whenever the function is re−issued,
    // : even if no leap seconds have been
    // : added.
    // :
    // : Latest leap second: 2016 December 31
    // : [NB: this is still true as of 2023-02-23]
    // :
    // :__________________________________________
    public static double Dat(int iy, int im, int id, double fd) {
        double deltat = 0;
        var status = S_Dat(iy, im, id, fd, ref deltat);
        return deltat; }
    [DllImport(DllFilename, EntryPoint = "iauDat", CallingConvention = CallingConvention.Cdecl)]
    static extern int S_Dat(int iy, int im, int id, double fd, ref double deltat);
    
    
    // /// <summary>Sofa.Dtdb(): An approximation to TDB−TT, the difference between barycentric
    // ///                       dynamical time and terrestrial time, for an observer on the Earth.</summary>
    // /// <param name="date1">Part 1 of 2-part date, TDB</param>
    // /// <param name="date2">Part 2 of 2-part date, TDB</param>
    // /// <param name="ut">Universal time (UT1, fraction of one day)</param>
    // /// <param name="elong">Longitude (radians, east positive)</param>
    // /// <param name="u">Distance from Earth spin axis (km)</param>
    // /// <param name="v">Distance north of equatorial plane (km) </param>
    // /// <returns>double: TAI minus UTC (seconds)</returns>
    // //-----------------------------------------------------------------
    // //   TAI                    <− physically realized
    // //    |
    // //   offset                 <− observed (nominally +32.184s)
    // //    |
    // //   TT                     <− terrestrial time
    // //    |
    // //   rate adjustment (L_G)  <− definition of TT
    // //    |
    // //   TCG                    <− time scale for GCRS
    // //    |
    // //   "periodic" terms       <− iauDtdb is an implementation
    // //    |
    // //   rate adjustment (L_C)  <− function of solar−system ephemeris
    // //    |
    // //   TCB                    <− time scale for BCRS
    // //    |
    // //   rate adjustment (−L_B) <− definition of TDB
    // //    |
    // //   TDB                    <− TCB scaled to track TT
    // //    |
    // //   "periodic" terms       <− −iauDtdb is an approximation
    // //    |
    // //   TT                     <− terrestrial time
    // //-----------------------------------------------------------------
    // // Adopted values for the various constants can be found in the
    // //     IERS Conventions (McCarthy & Petit 2003).
    // //-----------------------------------------------------------------
    // public static double Dtdb(double date1, double date2, double ut, double elong, double u, double v) {
    //     var diff = S_Dtdb(date1, date2, ut, elong, u, v);
    //     return diff; }
    // [DllImport(DllFilename, EntryPoint = "iauDtdb", CallingConvention = CallingConvention.Cdecl)]
    // static extern double S_Dtdb(double date1, double date2, double ut, double elong, double u, double v);
    
    
    /// <summary>Sofa.Dtf2d(): Encode date and time fields into 2−part Julian Date (or in the case
    ///                        of UTC a quasi−JD form that includes special provision for
    ///                        leap seconds).</summary>
    /// <param name="scale">Timescale ID, typically "UTC"</param>
    /// <param name="iy">Year in Gregorian calendar</param>
    /// <param name="im">Month in Gregorian calendar</param>
    /// <param name="id">Day in Gregorian calendar</param>
    /// <param name="ihr">Hour, UTC</param>
    /// <param name="imn">Minute, UTC</param>
    /// <param name="sec">Seconds, UTC</param>
    /// <returns>2-Tuple of doubles: d1, d2, a 2−part Julian Date</returns>
    public static Tuple<double, double> Dtf2d(string scale, int iy, int im, int id, 
        int ihr, int imn, double sec) {
        double d1 = 0, d2 = 0;
        var status = S_Dtf2d(scale, iy, im, id, ihr, imn, sec, ref d1, ref d2);
        return Tuple.Create(d1, d2); }
    [DllImport(DllFilename, EntryPoint = "iauDtf2d", CallingConvention = CallingConvention.Cdecl)]
    static extern int S_Dtf2d(string scale, int iy, int im, int id, int ihr, int imn, double sec, 
        ref double d1, ref double d2);
    
    
    /// <summary>Sofa.Taitt(): Time scale transformation: International Atomic Time, TAI,
    ///                        to Terrestrial Time, TT.</summary>
    /// <param name="tai1">Part 1 of 2-part Julian Date, TAI.</param>
    /// <param name="tai2">Part 2 of 2-part Julian Date, TAI.</param>
    /// <returns>2-Tuple of doubles: tt1: Part 1 of 2-part Julian Date, TT
    ///                              tt2: Part 2 of 2-part Julian Date, TT.</returns>
    public static Tuple<double, double> Taitt(double tai1, double tai2) {
        double tt1 = 0, tt2 = 0;
        var status = S_Taitt(tai1, tai2, ref tt1, ref tt2);
        return Tuple.Create(tt1, tt2); }
    [DllImport(DllFilename, EntryPoint = "iauTaitt", CallingConvention = CallingConvention.Cdecl)]
    static extern int S_Taitt(double tai1, double tai2, ref double tt1, ref double tt2);

    
    /// <summary>Sofa.Taiut1(): Time scale transformation: International Atomic Time, TAI,
    ///                         to Universal Time, UT1.</summary>
    /// <param name="tai1">Part 1 of 2-part Julian Date, TAI.</param>
    /// <param name="tai2">Part 2 of 2-part Julian Date, TAI.</param>
    /// <param name="dta">UT1-TAI (seconds).</param>
    /// <returns>2-Tuple of doubles: ut11: Part 1 of 2-part Julian Date, UT1
    ///                              ut12: Part 2 of 2-part Julian Date, UT1.</returns>
    public static Tuple<double, double> Taiut1(double tai1, double tai2, double dta) {
        double ut11 = 0, ut12 = 0;
        var status = S_Taiut1(tai1, tai2, dta, ref ut11, ref ut12);
        return Tuple.Create(ut11, ut12); }
    [DllImport(DllFilename, EntryPoint = "iauTaiut1", CallingConvention = CallingConvention.Cdecl)]
    static extern int S_Taiut1(double tai1, double tai2, double dta, ref double ut11, ref double ut12);

    
    /// <summary>Sofa.Taiutc(): Time scale transformation: International Atomic Time, TAI,
    ///                         to Coordinated Universal Time, UTC.</summary>
    /// <param name="tai1">Part 1 of 2-part Julian Date, TAI.</param>
    /// <param name="tai2">Part 2 of 2-part Julian Date, TAI.</param>
    /// <returns>2-Tuple of doubles: utc1: Part 1 of 2-part Julian Date, UTC
    ///                              utc2: Part 2 of 2-part Julian Date, UTC.</returns>
    public static Tuple<double, double> Taiutc(double tai1, double tai2) {
        double utc1 = 0, utc2 = 0;
        var status = S_Taiutc(tai1, tai2, ref utc1, ref utc2);
        return Tuple.Create(utc1, utc2); }
    [DllImport(DllFilename, EntryPoint = "iauTaiutc", CallingConvention = CallingConvention.Cdecl)]
    static extern int S_Taiutc(double tai1, double tai2, ref double utc1, ref double utc2);

    
    /// <summary>Sofa.Tcbtdb(): Time scale transformation: Barycentric Coordinate Time, TCB,
    ///                         to Barycentric Dynamical Time, TDB (~= Teph for JPL ephemerides).</summary>
    /// <param name="tcb1">Part 1 of 2-part Julian Date, TCB.</param>
    /// <param name="tcb2">Part 2 of 2-part Julian Date, TCB.</param>
    /// <returns>2-Tuple of doubles: tdb1: Part 1 of 2-part Julian Date, TDB
    ///                              tdb2: Part 2 of 2-part Julian Date, TDB.</returns>
    public static Tuple<double, double> Tcbtdb(double tcb1, double tcb2) {
        double tdb1 = 0, tdb2 = 0;
        var status = S_Tcbtdb(tcb1, tcb2, ref tdb1, ref tdb2);
        return Tuple.Create(tdb1, tdb2); }
    [DllImport(DllFilename, EntryPoint = "iauTcbtdb", CallingConvention = CallingConvention.Cdecl)]
    static extern int S_Tcbtdb(double tcb1, double tcb2, ref double tdb1, ref double tdb2);
    
    
    /// <summary>Sofa.Tcgtt(): Time scale transformation: Geocentric Coordinate Time, TCG,
    ///                        to Terrestrial Time, TT.</summary>
    /// <param name="tcg1">Part 1 of 2-part Julian Date, TCG.</param>
    /// <param name="tcg2">Part 2 of 2-part Julian Date, TCG.</param>
    /// <returns>2-Tuple of doubles: tt1: Part 1 of 2-part Julian Date, TT
    ///                              tt2: Part 2 of 2-part Julian Date, TT.</returns>
    public static Tuple<double, double> Tcgtt(double tcg1, double tcg2) {
        double tt1 = 0, tt2 = 0;
        var status = S_Tcgtt(tcg1, tcg2, ref tt1, ref tt2);
        return Tuple.Create(tt1, tt2); }
    [DllImport(DllFilename, EntryPoint = "iauTcgtt", CallingConvention = CallingConvention.Cdecl)]
    static extern int S_Tcgtt(double tcg1, double tcg2, ref double tt1, ref double tt2);

    
    /// <summary>Sofa.Tdbtcb(): Time scale transformation: Barycentric Dynamical Time, TDB,
    ///                         to Barycentric Coordinate Time, TCB.</summary>
    /// <param name="tdb1">Part 1 of 2-part Julian Date, TDB.</param>
    /// <param name="tdb2">Part 2 of 2-part Julian Date, TDB.</param>
    /// <returns>2-Tuple of doubles: tcb1: Part 1 of 2-part Julian Date, TCB
    ///                              tcb2: Part 2 of 2-part Julian Date, TCB.</returns>
    public static Tuple<double, double> Tdbtcb(double tdb1, double tdb2) {
        double tcb1 = 0, tcb2 = 0;
        var status = S_Tdbtcb(tdb1, tdb2, ref tcb1, ref tcb2);
        return Tuple.Create(tcb1, tcb2); }
    [DllImport(DllFilename, EntryPoint = "iauTdbtcb", CallingConvention = CallingConvention.Cdecl)]
    static extern int S_Tdbtcb(double tdb1, double tdb2, ref double tcb1, ref double tcb2);

            
    /// <summary>Sofa.Tttai(): Time scale transformation: Terrestrial Time, TT,
    ///                        to International Atomic Time, TAI.</summary>
    /// <param name="tt1">Part 1 of 2-part Julian Date, TT.</param>
    /// <param name="tt2">Part 2 of 2-part Julian Date, TT.</param>
    /// <returns>2-Tuple of doubles: tai1: Part 1 of 2-part Julian Date, TAI
    ///                              tai2: Part 2 of 2-part Julian Date, TAI.</returns>
    public static Tuple<double, double> Tttai(double tt1, double tt2) {
        double tai1 = 0, tai2 = 0;
        var status = S_Tttai(tt1, tt2, ref tai1, ref tai2);
        return Tuple.Create(tai1, tai2); }
    [DllImport(DllFilename, EntryPoint = "iauTttai", CallingConvention = CallingConvention.Cdecl)]
    static extern int S_Tttai(double tt1, double tt2, ref double tai1, ref double tai2);

            
    /// <summary>Sofa.Tttcg(): Time scale transformation: Terrestrial Time, TT,
    ///                        to Geocentric Coordinate Time, TCG.</summary>
    /// <param name="tt1">Part 1 of 2-part Julian Date, TT.</param>
    /// <param name="tt2">Part 2 of 2-part Julian Date, TT.</param>
    /// <returns>2-Tuple of doubles: tcg1: Part 1 of 2-part Julian Date, TCG
    ///                              tcg2: Part 2 of 2-part Julian Date, TCG.</returns>
    public static Tuple<double, double> Tttcg(double tt1, double tt2) {
        double tcg1 = 0, tcg2 = 0;
        var status = S_Tttcg(tt1, tt2, ref tcg1, ref tcg2);
        return Tuple.Create(tcg1, tcg2); }
    [DllImport(DllFilename, EntryPoint = "iauTttcg", CallingConvention = CallingConvention.Cdecl)]
    static extern int S_Tttcg(double tt1, double tt2, ref double tcg1, ref double tcg2);

    
    /// <summary>Sofa.Ttut1(): Time scale transformation: Terrestrial Time, TT,
    ///                        to Universal Time, UT1.</summary>
    /// <param name="tt1">Part 1 of 2-part Julian Date, TT.</param>
    /// <param name="tt2">Part 2 of 2-part Julian Date, TT.</param>
    /// <param name="dt">TT-UT1 (seconds).</param>
    /// <returns>2-Tuple of doubles: tut11: Part 1 of 2-part Julian Date, UT1
    ///                              tut12: Part 2 of 2-part Julian Date, UT1.</returns>
    public static Tuple<double, double> Ttut1(double tt1, double tt2, double dt) {
        double tut11 = 0, tut12 = 0;
        var status = S_Ttut1(tt1, tt2, dt, ref tut11, ref tut12);
        return Tuple.Create(tut11, tut12); }
    [DllImport(DllFilename, EntryPoint = "iauTtut1", CallingConvention = CallingConvention.Cdecl)]
    static extern int S_Ttut1(double tt1, double tt2, double dt, ref double tut11, ref double tut12);

    
    /// <summary>Sofa.Ut1tai(): TTime scale transformation: Universal Time, UT1,
    ///                         to International Atomic Time, TAI.</summary>
    /// <param name="ut11">Part 1 of 2-part Julian Date, UT1.</param>
    /// <param name="ut12">Part 2 of 2-part Julian Date, UT1.</param>
    /// <param name="dta">UT1-TAI, available IERS tabulations (seconds).</param>
    /// <returns>2-Tuple of doubles: tai1: Part 1 of 2-part Julian Date, TAI
    ///                              tai2: Part 2 of 2-part Julian Date, TAI.</returns>
    public static Tuple<double, double> Ut1tai(double ut11, double ut12, double dta) {
        double tai1 = 0, tai2 = 0;
        var status = S_Ut1tai(ut11, ut12, ref tai1, ref tai2);
        return Tuple.Create(tai1, tai2); }
    [DllImport(DllFilename, EntryPoint = "iauUt1tai", CallingConvention = CallingConvention.Cdecl)]
    static extern int S_Ut1tai(double ut11, double ut12, ref double tai1, ref double tai2);

        
    /// <summary>Sofa.Ut1tt(): Time scale transformation: Universal Time, UT1,
    ///                        to Terrestrial Time, TT.</summary>
    /// <param name="ut11">Part 1 of 2-part Julian Date, UT1.</param>
    /// <param name="ut12">Part 2 of 2-part Julian Date, UT1.</param>
    /// <param name="dt">TT-UT1, "classical Delta T" (seconds).</param>
    /// <returns>2-Tuple of doubles: tt1: Part 1 of 2-part Julian Date, TT
    ///                              tt2: Part 2 of 2-part Julian Date, TT.</returns>
    public static Tuple<double, double> Ut1tt(double ut11, double ut12, double dt) {
        double tt1 = 0, tt2 = 0;
        var status = S_Ut1tt(ut11, ut12, dt, ref tt1, ref tt2);
        return Tuple.Create(tt1, tt2); }
    [DllImport(DllFilename, EntryPoint = "iauUt1tt", CallingConvention = CallingConvention.Cdecl)]
    static extern int S_Ut1tt(double ut11, double ut12, double dt, ref double tt1, ref double tt2);

    
    /// <summary>Sofa.Ut1utc(): Time scale transformation: Universal Time, UT1,
    ///                         to Coordinated Universal Time, UTC.</summary>
    /// <param name="ut11">Part 1 of 2-part Julian Date, UT1.</param>
    /// <param name="ut12">Part 2 of 2-part Julian Date, UT1.</param>
    /// <param name="dut1">UT1-UTC, "Delta UT1", available IERS tabulations (seconds).</param>
    /// <returns>2-Tuple of doubles: utc1: Part 1 of 2-part Julian Date, UTC
    ///                              utc2: Part 2 of 2-part Julian Date, UTC.</returns>
    public static Tuple<double, double> Ut1utc(double ut11, double ut12, double dut1) {
        double utc1 = 0, utc2 = 0;
        var status = S_Ut1utc(ut11, ut12, dut1, ref utc1, ref utc2);
        return Tuple.Create(utc1, utc2); }
    [DllImport(DllFilename, EntryPoint = "iauUt1utc", CallingConvention = CallingConvention.Cdecl)]
    static extern int S_Ut1utc(double ut11, double ut12, double dut1, ref double utc1, ref double utc2);

    
    /// <summary>Sofa.Utctai(): Time scale transformation: Coordinated Universal Time, UTC,
    ///                         to International Atomic Time, TAI.</summary>
    /// <param name="utc1">Part 1 of 2-part Julian Date, UTC.</param>
    /// <param name="utc2">Part 2 of 2-part Julian Date, UTC.</param>
    /// <returns>2-Tuple of doubles: tai1: Part 1 of 2-part Julian Date, TAI
    ///                              tai2: Part 2 of 2-part Julian Date, TAI.</returns>
    public static Tuple<double, double> Utctai(double utc1, double utc2) {
        double tai1 = 0, tai2 = 0;
        var status = S_Utctai(utc1, utc2, ref tai1, ref tai2);
        return Tuple.Create(tai1, tai2); }
    [DllImport(DllFilename, EntryPoint = "iauUtctai", CallingConvention = CallingConvention.Cdecl)]
    static extern int S_Utctai(double utc1, double utc2, ref double tai1, ref double tai2);

    
    /// <summary>Sofa.Utcut1(): Time scale transformation: Coordinated Universal Time, UTC,
    ///                         to Universal Time, UT1.</summary>
    /// <param name="utc1">Part 1 of 2-part Julian Date, UTC.</param>
    /// <param name="utc2">Part 2 of 2-part Julian Date, UTC.</param>
    /// <param name="dut1">UT1-UTC, "Delta UT1", available IERS tabulations (seconds).</param>
    /// <returns>2-Tuple of doubles: ut11: Part 1 of 2-part Julian Date, UT1
    ///                              ut12: Part 2 of 2-part Julian Date, UT1.</returns>
    public static Tuple<double, double> Utcut1(double utc1, double utc2, double dut1) {
        double ut11 = 0, ut12 = 0;
        var status = S_Utcut1(utc1, utc2, dut1, ref ut11, ref ut12);
        return Tuple.Create(ut11, ut12); }
    [DllImport(DllFilename, EntryPoint = "iauUtcut1", CallingConvention = CallingConvention.Cdecl)]
    static extern int S_Utcut1(double utc1, double utc2, double dut1, ref double ut11, ref double ut12);
#endregion    
    
#region "Astronomy/HorizonEquatorial"
    /*######### Astronomy/HorizonEquatorial ##############################################################*/
    // No need to expose as public *for now* (reproduce these as needed, in C# within AstroLib.Geometry):
    //      iauHd2pa  Parallactic angle for a given hour angle and declination.
    // .....................................................................................................
    
    
    /// <summary>Sofa.Ae2hd(): Horizon to equatorial coordinates: transform azimuth and altitude to
    ///                        hour angle and declination.</summary>
    /// <param name="az">Azimuth (radians).</param>
    /// <param name="el">Altitude (elevation, radians).</param>
    /// <param name="phi">Site latitude (radians).</param>
    /// <returns>2-tuple of doubles:
    ///              ha: Local hour angle (radians)
    ///              dec: Site declination (radians)</returns>.
    public static Tuple<double, double> Ae2hd(double az, double el, double phi) {
        double ha = 0, dec = 0;
        S_Ae2hd(az, el, phi, ref ha, ref dec);
        return Tuple.Create(ha, dec); }
    [DllImport(DllFilename, EntryPoint = "iauAe2hd", CallingConvention = CallingConvention.Cdecl)]
    static extern int S_Ae2hd(double az, double el, double phi, ref double ha, ref double dec); 
 
    
    /// <summary>Sofa.Hd2ae(): Equatorial to horizon coordinates: transform hour angle and
    ///                        declination to azimuth and altitude.</summary>
    /// <param name="ha">Local hour angle (radians).</param>
    /// <param name="dec">Site declination (radians).</param>
    /// <param name="phi">Site latitude (radians).</param>
    /// <returns>2-tuple of doubles:
    ///              az: Azimuth (radians)
    ///              el: Altitude (elevation, radians)</returns>.
    public static Tuple<double, double> Hd2ae(double ha, double dec, double phi) {
        double az = 0, el = 0;
        S_Hd2ae(ha, dec, phi, ref az, ref el);
        return Tuple.Create(az, el); }
    [DllImport(DllFilename, EntryPoint = "iauHd2ae", CallingConvention = CallingConvention.Cdecl)]
    static extern int S_Hd2ae(double ha, double dec, double phi, ref double az, ref double el);

#endregion

#region "Astronomy/Gnomonic"
    /*######### Astronomy/Gnomonic ######################################################################*/
    // No need to expose as public *for now* (reproduce these as needed, in C# within AstroLib.Geometry):
    //      iauTpors  In the tangent plane projection, given the rectangular coordinates of a star and its
    //                    spherical coordinates, determine the spherical coordinates of the tangent point.
    //      iauTporv  In the tangent plane projection, given the rectangular coordinates of a star and its
    //                    direction cosines, determine the direction cosines of the tangent point.
    //      iauTpstv  In the tangent plane projection, given the star’s rectangular coordinates and the
    //                    direction cosines of the tangent point, solve for the direction cosines
    //                    of the star.
    //      iauTpxev  In the tangent plane projection, given celestial direction cosines for a star and
    //                    the tangent point, solve for the star’s rectangular coordinates
    //                    in the tangent plane.
    // .....................................................................................................

    
    /// <summary>Sofa.Tpsts(): In the tangent plane projection, given the star’s rectangular coordinates
    ///                        and the spherical coordinates of the tangent point, solve for the
    ///                        spherical coordinates of the star.</summary>
    /// <param name="xi">Rectangular coordinates of star image (radians at tangent point).</param>
    /// <param name="eta">Rectangular coordinates of star image, due North (radians @ tangent pt).</param>
    /// <param name="a0">Tangent point's spherical coordinates (radians).</param>
    /// <param name="b0">Tangent point's spherical coordinates (radians).</param>
    /// <returns>2-tuple of doubles:
    ///              a: Star's spherical coordinates (radians).
    ///              b: Star's spherical coordinates (radians).</returns>.
    public static Tuple<double, double> Tpsts(double xi, double eta, double a0, double b0) {
        double a = 0, b = 0;
        S_Tpsts(xi, eta, a0, b0, ref a, ref b);
        return Tuple.Create(a, b); }
    [DllImport(DllFilename, EntryPoint = "iauTpsts", CallingConvention = CallingConvention.Cdecl)]
    static extern int S_Tpsts(double xi, double eta, double a0, double b0, ref double a, ref double b);

    
    /// <summary>Sofa.Tpxes(): In the tangent plane projection, given celestial spherical coordinates for
    ///                        a star and the tangent point, solve for the star’s rectangular coordinates
    ///                        in the tangent plane.</summary>
    /// <param name="a">Star's spherical coordinates (radians).</param>
    /// <param name="b">Star's spherical coordinates (radians).</param>
    /// <param name="a0">Tangent point's spherical coordinates (radians).</param>
    /// <param name="b0">Tangent point's spherical coordinates (radians).</param>
    /// <returns>2-tuple of doubles:
    ///                       xi:  Rectangular coordinates of star image (radians at tangent point).
    ///                       eta: Rectangular coordinates of star image, due North (radians at tangent pt).
    ///                       </returns>.
    public static Tuple<double, double> Tpxes(double a, double b, double a0, double b0) {
        double xi = 0, eta = 0;
        var status = S_Tpxes(a, b, a0, b0, ref xi, ref eta);
        return Tuple.Create(xi, eta); }
    [DllImport(DllFilename, EntryPoint = "iauTpxes", CallingConvention = CallingConvention.Cdecl)]
    static extern int S_Tpxes(double a, double b, double a0, double b0, ref double xi, ref double eta);

    
    
    
#endregion

#region "VectorMatrix (all)"
    /*######### VectorMatrix/AngleOps ####################################################################*/
    // No need to expose as public *for now* (reproduce these as needed, in C# within AstroLib.Geometry):
    //      iauAf2a: Convert degrees, arcminutes, arcseconds to radians.
    //      iauAnp:  Normalize angle into the range 0 <= a < 2pi.
    //      iauAnpm: Normalize angle into the range -pi <= a < +pi.
    //      iauD2tf: Decompose days to hours, minutes, seconds, fraction.
    //      iauTf2d   Convert hours, minutes, seconds to days.
    // .....................................................................................................
    
    
    /// <summary> iauA2af: Decompose radians into degrees, arcminutes, arcseconds, fraction.</summary>
    /// <param name="ndp">Resolution, as negative 10-exponent of last decimal place,
    /// e.g., -2 for rounding to hundreds, +1 for rounding to 1/10, +6 for rounding to millionths.</param>
    /// <param name="angle">Input angle (radians). May exceed 2*pi.
    /// NB: user is responsible for handling when near 360 degrees.</param>
    /// <returns>2-Tuple:
    ///     sign (string of length one): '+' or '-'.
    ///     idmsf (int array of length 4):
    ///         degrees, arcminutes, arcseconds, fraction in the requested resolution.</returns>
    // TODO: In testing: determine what idmsf[3] (fraction) is supposed to represent for various values of ndp.
    public static Tuple<string, int[]> A2af(int ndp, double angle) {
        string sign = " ";
        var idmsf = new int[4] {0, 0, 0, 0}; 
        S_A2af(ndp, angle, sign, idmsf);
        var result = Tuple.Create(sign, idmsf);
        return result;
    }
    [DllImport(DllFilename, EntryPoint = "iauA2af", CallingConvention = CallingConvention.Cdecl)]
    static extern void S_A2af(int ndp, double angle, string sign, int[] idmsf);
    
    
    /// <summary> iauA2tf: Decompose radians into hours, minutes, seconds, fraction.</summary>
    /// <param name="ndp">Resolution, as negative 10-exponent of last decimal place,
    /// e.g., -2 for rounding to hundreds, +1 for rounding to 1/10, +6 for rounding to millionths.</param>
    /// <param name="angle">Input angle (radians). May exceed 2*pi.
    /// NB: user is responsible for handling when near 360 degrees.</param>
    /// <returns>2-Tuple:
    ///     sign (string of length one): '+' or '-'.
    ///     ihmsf (int array of length 4):
    ///         hours, minutes, seconds, fraction in the requested resolution.</returns>
    // TODO: Determine what ihmsf[3] (fraction) is supposed to represent for various values of ndp.
    public static Tuple<string, int[]> A2tf(int ndp, double angle) {
        string sign = " ";
        var ihmsf = new int[4] {0, 0, 0, 0}; 
        S_A2tf(ndp, angle, sign, ihmsf);
        var result = Tuple.Create(sign, ihmsf);
        return result;
    }
    [DllImport(DllFilename, EntryPoint = "iauA2tf", CallingConvention = CallingConvention.Cdecl)]
    static extern void S_A2tf(int ndp, double angle, string sign, int[] ihmsf);


    /// <summary> Sofa.Tf2a(): Convert hours, minutes, seconds to angle in radians.</summary>
    /// <param name="s">sign (char), '-' means negative, any other means positive.</param>
    /// <param name="ihour">hours (int)</param>
    /// <param name="imin">minutes (int)</param>
    /// <param name="sec">seconds (double)</param>
    /// <returns>angle (radians)</returns>
    public static double Tf2a(char s, int ihour, int imin, double sec) {
        double rad = 0;
        var status = S_Tf2a(s, ihour, imin, sec, ref rad);
        return rad;
    }
    [DllImport(DllFilename, EntryPoint = "iauTf2a", CallingConvention = CallingConvention.Cdecl)]
    static extern int S_Tf2a(char s, int ihour, int imin, double sec, ref double rad);
    
    
    /// <summary> Sofa.Tf2d(): Convert hours, minutes, seconds to days. </summary>
    /// <param name="s">sign (char), '-' means negative, any other means positive.</param>
    /// <param name="ihour">hours (int)</param>
    /// <param name="imin">minutes (int)</param>
    /// <param name="sec">seconds (double)</param>
    /// <returns>days (double)</returns>
    public static double Tf2d(char s, int ihour, int imin, double sec) {
        double days = 0;
        var status = S_Tf2a(s, ihour, imin, sec, ref days);
        return days;
    }
    [DllImport(DllFilename, EntryPoint = "iauTf2d", CallingConvention = CallingConvention.Cdecl)]
    static extern int S_Tf2d(char s, int ihour, int imin, double sec, ref double days);
    
    
    /*######### VectorMatrix/BuildRotations ##############################################################*/
    // .....................................................................................................


    /// <summary> Sofa.Rx(): Rotate a rotation matrix about the x−axis.</summary>
    /// <param name="phi">Newly applied rotation angle (radians).</param>
    /// <param name="r">Rotation matrix before new rotation.</param>
    /// <returns>Rotation matrix after new rotation.</returns>
    public static double[,] Rx(double phi, double[,] r) {
        S_Rx(phi, r);
        return r;
    }
    [DllImport(DllFilename, EntryPoint = "iauRx", CallingConvention = CallingConvention.Cdecl)]
    static extern void S_Rx(double phi, [In, Out] double[,] r);
    
    
    /// <summary> Sofa.Ry(): Rotate a rotation matrix about the y−axis.</summary>
    /// <param name="theta">Newly applied rotation angle (radians).</param>
    /// <param name="r">Rotation matrix before new rotation.</param>
    /// <returns>Rotation matrix after new rotation.</returns>
    public static double[,] Ry(double theta, double[,] r) {
        S_Ry(theta, r);
        return r;
    }
    [DllImport(DllFilename, EntryPoint = "iauRy", CallingConvention = CallingConvention.Cdecl)]
    static extern void S_Ry(double theta, [In, Out] double[,] r);
    
    
    /// <summary> Sofa.Rz(): Rotate a rotation matrix about the z−axis.</summary>
    /// <param name="psi">Newly applied rotation angle (radians).</param>
    /// <param name="r">Rotation matrix before new rotation.</param>
    /// <returns>Rotation matrix after new rotation.</returns>
    public static double[,] Rz(double psi, double[,] r) {
        S_Rz(psi, r);
        return r;
    }
    [DllImport(DllFilename, EntryPoint = "iauRz", CallingConvention = CallingConvention.Cdecl)]
    static extern void S_Rz(double psi, [In, Out] double[,] r);

    
    /*######### VectorMatrix/CopyExtendExtract ###########################################################*/
    // No need to expose as public *for now* (reproduce these as needed, in C# within AstroLib.Geometry):
    //      iauCp:   Copy a p-vector.
    //      iauCpv:  Copy a position/velocity vector.
    //      iauCr:   Copy an r-matrix.
    //      iauP2pv   Extend a p−vector to a pv−vector by appending a zero velocity.
    //      iauPv2p   Discard velocity component of a pv−vector.
    // .....................................................................................................


    /*######### VectorMatrix/Initialization ##############################################################*/
    // No need to expose as public *for now* (reproduce these as needed, in C# within AstroLib.Geometry):
    //      iauIr:    Initialize an r−matrix to the identity matrix.    
    //      iauZp     Zero a p−vector.
    //      iauZpv    Zero a pv−vector.
    //      iauZr     Initialize an r−matrix to the null matrix.
    // .....................................................................................................
    
    
    /*######### VectorMatrix/MatrixOps ###################################################################*/
    // No need to expose as public *for now* (reproduce these as needed, in C# within AstroLib.Geometry):
    // .....................................................................................................


    /// <summary> Multiply two rotation matrices. NB: not commutative, Rxr(a, b) != Rxr(b, a).
    /// Same array may be used for any of the arguments.</summary>
    /// <param name="a">First rotation matrix.</param>
    /// <param name="b">Second rotation matrix.</param>
    /// <returns>a * b.</returns>
    public static double[,] Rxr(double[,] a, double[,] b) {
        var atb = new double[3, 3] {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
        S_Rxr(a, b, atb);
        return atb;
    }
    [DllImport(DllFilename, EntryPoint = "iauRxr", CallingConvention = CallingConvention.Cdecl)]
    static extern void S_Rxr(double[,] a, double[,] b, [In, Out] double[,] atb);
    
    
    /// <summary> Return transpose of a rotation matrix.</summary>
    /// <param name="r">Rotation matrix before transposition.</param>
    /// <returns>Transposed rotation matrix.</returns>
    public static double[,] Tr(double[,] r) {
        var rt = new double[3, 3] {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
        S_Tr(r, rt);
        return rt;
    }
    [DllImport(DllFilename, EntryPoint = "iauTr", CallingConvention = CallingConvention.Cdecl)]
    static extern void S_Tr(double[,] r, [In, Out] double[,] rt);
    
    
    /*######### VectorMatrix/MatrixVectorProducts ########################################################*/
    // No need to expose as public *for now* (reproduce these as needed, in C# within AstroLib.Geometry):
    //      *** iauRxp    Multiply a p−vector by an r−matrix.
    //      iauRxpv   Multiply a pv−vector by an r−matrix.
    //      iauTrxp   Multiply a p−vector by the transpose of an r−matrix.
    //      iauTrxpv  Multiply a pv−vector by the transpose of an r−matrix.
    // .....................................................................................................


    /// <summary>Sofa.Rxp(): Multiply a position vector by a rotation matrix.</summary>
    /// <param name="r">Rotation matrix.</param>
    /// <param name="p">Position vector.</param>
    /// <returns>Position vector after rotation, i.e., r * p.</returns>
    public static double[] Rxp(double[,] r, double[] p) {
        var rp = new double[3] {0, 0, 0};
        S_Rxp(r, p, rp);
        return rp;
    }
    [DllImport(DllFilename, EntryPoint = "iauRxp", CallingConvention = CallingConvention.Cdecl)]
    static extern void S_Rxp(double[,] r, double[] p, [In, Out] double[] rp);

    
    /*######### VectorMatrix/RotationVectors #############################################################*/
    // No need to expose as public *for now* (reproduce these as needed, in C# within AstroLib.Geometry):
    //      iauRm2v   Express an r−matrix (xyz) as an r−vector (per Euler axis).
    //      iauRv2m   Form the r−matrix corresponding to a given r−vector.
    // .....................................................................................................
    
    
    /*########## VectorMatrix/SeparationAndAngle #########################################################*/
    // No need to expose as public *for now* (reproduce these as needed, in C# within AstroLib.Geometry):
    //      iauPap    Position−angle from two p−vectors.
    //      iauPas    Position−angle from spherical coordinates.
    // .....................................................................................................
    

    /// <summary>Sofa.Sepp(): Angular separation between two p−vectors.</summary>
    /// <param name="a">First p−vector (not necessarily unit length)</param>
    /// <param name="b">Second p−vector (not necessarily unit length)</param>
    /// <returns>Angular separation (radians, always in range [0, pi])</returns>
    public static double Sepp(double[] a, double[] b) {
        return S_Sepp(a, b); }
    [DllImport(DllFilename, EntryPoint = "iauSepp", CallingConvention = CallingConvention.Cdecl)]
    static extern double S_Sepp(double[] a, double[] b);
    
    
    /// <summary>Sofa.Seps(): Angular separation between two sets of spherical coordinates.</summary>
    /// <param name="al">First longitude/RA (radians)</param>
    /// <param name="ap">First latitude/Declination (radians)</param>
    /// <param name="bl">Second longitude/RA (radians)</param>
    /// <param name="bp">Second latitude/Declination (radians)</param>
    /// <returns>Angular separation (radians)</returns>
    public static double Seps(double al, double ap, double bl, double bp) {
        return S_Seps(al, ap, bl, bp); }
    [DllImport(DllFilename, EntryPoint = "iauSeps", CallingConvention = CallingConvention.Cdecl)]
    static extern double S_Seps(double al, double ap, double bl, double bp);
    
    
    /*######### VectorMatrix/SphericalCartesian ##########################################################*/
    // No need to expose as public *for now* (reproduce these as needed, in C# within AstroLib.Geometry):
    //      iauS2c    Convert spherical coordinates to Cartesian (direction cosines).
    // .....................................................................................................
    
    
    /// <summary>Sofa.C2s(): Convert position vector to spherical coordinates.</summary>
    /// <param name="p">Position vector (double[3]) of any magnitude, only its direction is used.</param>
    /// <returns>2-Tuple:
    ///     theta (double): longitude angle (radians)
    ///     phi (double): latitude angle (radians)</returns>
    public static Tuple<double, double> C2s(double[] p) {
        double theta = 0, phi = 0;
        S_C2s(p, ref theta, ref phi);
        var result = Tuple.Create(theta, phi);
        return result;
    }    
    [DllImport(DllFilename, EntryPoint = "iauC2s", CallingConvention = CallingConvention.Cdecl)]
    static extern double S_C2s(double[] p, ref double theta, ref double phi);
    
    
    /// <summary>Sofa.P2s(): Convert position vector to spherical polar coordinates.</summary>
    /// <param name="p">Position vector (double[3]) of any magnitude, only its direction is used.</param>
    /// <returns>3-Tuple:
    ///     theta (double): longitude angle (radians)
    ///     phi (double): latitude angle (radians)
    ///     r (double): radial distance</returns>
    public static Tuple<double, double, double> P2s(double[] p) {
        double theta = 0, phi = 0, r = 0;
        S_P2s(p, ref theta, ref phi, ref r);
        var result = Tuple.Create(theta, phi, r);
        return result;
    }    
    [DllImport(DllFilename, EntryPoint = "iauP2s", CallingConvention = CallingConvention.Cdecl)]
    static extern double S_P2s(double[] p, ref double theta, ref double phi, ref double r);

    
    /// <summary>Sofa.Pv2s(): Convert position/velocity from Cartesian to spherical coordinates
    /// The inverse of function Sofa.S2pv().</summary>
    /// <param name="pv">Position/velocity vectors as double[2,3]</param>
    /// <returns>6-Tuple:
    ///     theta (double): longitude angle (radians)
    ///     phi (double): latitude angle (radians)
    ///     r (double): radial distance
    ///     td (double): rate of change of theta
    ///     pd (double): rate of change of phi
    ///     rd (double): rate of change of radius</returns>
    public static Tuple<double, double, double, double, double, double> Pv2s(double[,] pv) {
        double theta = 0, phi = 0, r = 0, td = 0, pd = 0, rd = 0;
        S_Pv2s(pv, ref theta, ref phi, ref r, ref td, ref pd, ref rd);
        var result = Tuple.Create(theta, phi, r, td, pd, rd);
        return result;
    }    
    [DllImport(DllFilename, EntryPoint = "iauPv2s", CallingConvention = CallingConvention.Cdecl)]
    static extern double S_Pv2s(double[,] pv, ref double theta, ref double phi, ref double r,
        ref double td, ref double pd, ref double rd);


    /// <summary> Sofa.S2p(): Convert spherical polar coordinates to p−vector.</summary>
    /// <param name="theta">longitude angle (radians)</param>
    /// <param name="phi">latitude angle (radians)</param>
    /// <param name="r">radial distance</param>
    /// <returns>Cartesian coordinates.</returns>
    public static double[] S2p(double theta, double phi, double r) {
        var p = new double[3] {0, 0, 0};
        S_S2p(theta, phi, r, p);
        return p;
    }
    [DllImport(DllFilename, EntryPoint = "iauS2p", CallingConvention = CallingConvention.Cdecl)]
    static extern void S_S2p(double theta, double phi, double r, double[] p);
    

    /// <summary> Sofa.S2pv(): Convert position/velocity from spherical to Cartesian coordinates.
    /// The inverse of function Sofa.Pv2s().</summary>
    /// <param name="theta">longitude angle (radians)</param>
    /// <param name="phi">latitude angle (radians)</param>
    /// <param name="r">radial distance</param>
    /// <param name="td">rate of change of theta</param>
    /// <param name="pd">rate of change of phi</param>
    /// <param name="rd">rate of change of radius</param>
    /// <returns>Position/velicity vectors as double[2,3].</returns>
    public static double[,] S2pv(double theta, double phi, double r, double td, double pd, double rd) {
        var pv = new double[2, 3] {{0, 0, 0}, {0, 0, 0}};
        S_S2pv(theta, phi, r, td, pd, rd, pv);
        return pv;
    }
    [DllImport(DllFilename, EntryPoint = "iauS2pv", CallingConvention = CallingConvention.Cdecl)]
    static extern void S_S2pv(double theta, double phi, double r,
        double td, double pd, double rd, [In, Out] double[,] pv);

    
    /*######### VectorMatrix/VectorOps ###################################################################*/
    // No need to expose as public *for now* (reproduce these as needed, in C# within AstroLib.Geometry):
    //      iauPdp    p−vector inner (=scalar=dot) product.
    //      iauPm     Modulus of p−vector.
    //      iauPmp    P−vector subtraction.
    //      iauPn     Convert a p−vector into modulus and unit vector.
    //      iauPpp    P−vector addition.
    //      iauPpsp   P−vector plus scaled p−vector.
    //      iauPvdpv  Inner (=scalar=dot) product of two pv−vectors.
    //      iauPvm    Modulus of pv−vector.
    //      iauPvmpv  Subtract one pv−vector from another.
    //      iauPvppv  Add one pv−vector to another.
    //      iauPvu    Update (for new datetime) a pv−vector.
    //      iauPvup   Update (for new datetime) a pv−vector, discarding the velocity component.
    //      iauPvxpv  Outer (=vector=cross) product of two pv−vectors.
    //      iauPxp    p−vector outer (=vector=cross) product.
    //      iauS2xpv  Multiply a pv−vector by two scalars.
    //      iauSxp    Multiply a p−vector by a scalar.
    //      iauSxpv   Multiply a pv−vector by a scalar.
    // .....................................................................................................
#endregion
}