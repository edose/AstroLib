using System.Runtime.InteropServices;
using AstroLib.Core.Utils;

// ReSharper disable MemberCanBePrivate.Global
// ReSharper disable FieldCanBeMadeReadOnly.Global
// ReSharper disable StringLiteralTypo
// ReSharper disable CommentTypo
// ReSharper disable IdentifierTypo
// ReSharper disable InconsistentNaming

namespace AstroLib.Core.SOFA;

/*########################################################################################################*/
//
// This (contents of AstroLib.Core.SOFA) is not a source distribution or derived work of the
// IAU SOFA library. This is, rather, a set of our own methods that call SOFA functions.
// No SOFA source code is included here, though some SOFA documentation excerpts and SOFA-consistent
// parameter and struct naming is followed, to clarify use of AstroLib itself.
//
/*########################################################################################################*/
//
// Due diligence statement from the principal author of this class and its resulting assembly (.dll file):
// This work (AstroLib.Core.SOFA):
// (i) calls routines and computations that I derived by calling on but in no way reproducing software
//       provided by SOFA under license;
// (ii) does not itself constitute software provided by and/or endorsed by SOFA; and
// (iii) calls *compiled* C-language functions from C-language source code provided by SOFA,
//       but does not distribute that original source code.
//       We offer instead our own .NET-based API to access SOFA functions.
//       This new API is especially suitable for calling from another user's own C#-language code.
//
//       Eric Dose, Albuquerque, New Mexico USA 2023
//
/*########################################################################################################*/
//  Correspondence concerning SOFA software proper should be addressed as follows:
//
//      By email:  sofa@ukho.gov.uk
//      By post:   IAU SOFA Center
//                 HM Nautical Almanac Office
//                 UK Hydrographic Office
//                 Admiralty Way, Taunton
//                 Somerset, TA1 2DN
//                 United Kingdom
/*########################################################################################################*/

// Code regions in the following (e.g., "Astronomy/Calendars") group our AstroLib functions just as 
// the SOFA website (http://iausofa.org/2021_0512_C/Astronomy.html) groups the SOFA functions we call.

/// <summary>Class Sofa: Low-level .NET representations of many, but not all, IAU SOFA functions.
/// Contains functions from both the SOFA Astronomy and SOFA Vector/Matrix libraries. </summary>
public static partial class Sofa {
    private const string DllDirectory = "C:/DevCS/SOFA_dll/x64/Release";
    private const string DllFilename = "SOFA_dll.dll";

    static Sofa() {
        var dllFullpath = Path.Combine(DllDirectory, DllFilename);
        DllManager.LoadDllFromFullpath(dllFullpath);
    }



    
    // *****************************************************************************************************
    // Code regions in the following (e.g., "Astronomy: Calendars") group our AstroLib.Core.Sofa functions
    // just as the SOFA website (http://iausofa.org/2021_0512_C/Astronomy.html) groups them.
    // *****************************************************************************************************

    #region "Astronomy: Astrometry"
    /*######### Astronomy: Astrometry ####################################################################*/
    
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
    /// <param name="pnat">Natural direction to the source (unit vector)</param>
    /// <param name="v">Observer barycentric velocity in units of c</param>
    /// <param name="s">Distance between the Sun and the observer (AU)</param>
    /// <param name="bm1">Sqrt(1 - |v|^2), the reciprocal of Lorenz factor</param>
    /// <returns>Proper direction to the source (unit vector)</returns>
    public static double[] Ab(double[] pnat, double[] v, double s, double bm1) {
        var ppr = new double[3] {0, 0, 0};
        S_Ab(pnat, v, s, bm1, ppr);
        return ppr; }
    [DllImport(DllFilename, EntryPoint = "iauAb", CallingConvention = CallingConvention.Cdecl)]
    static extern void S_Ab(double[] pnat, double[] v, double s, double bm1, [In, Out] double[] ppr);
    
    // Apcg():  NOT IMPLEMENTED: earth ephemeris must be supplied by caller, so prefer Apcg13().

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
        return astrom; }
    [DllImport(DllFilename, EntryPoint = "iauApcg13", CallingConvention = CallingConvention.Cdecl)]
    static extern void S_Apcg13(double date1, double date2, ref iauASTROM astrom);

    // Apci():  NOT IMPLEMENTED: earth ephemeris must be supplied by caller, so prefer Apci13(). 

    /// <summary> Sofa.Apci13(): For a TERRESTRIAL OBSERVER, prepare star−independent astrometry parameters
    /// for transformations between ICRS and GCRS coordinates. The caller supplies the date,
    /// and SOFA models are used to predict the Earth ephemeris. The parameters produced by this function
    /// are required in the parallax, light deflection and aberration parts of the astrometric
    /// transformation chain.</summary>
    /// <param name="date1">TDB as Julian date, part 1.
    /// TT may be used without significant loss of accuracy.</param>
    /// <param name="date2">TDB as Julian date, part 2.
    /// TT may be used without significant loss of accuracy.</param>
    /// <returns> 2-Tuple:
    ///     astrom: Star−independent astrometry parameters, as a iauASTROM struct.
    ///             Unchanged struct elements: along, xp1, yp1, sphi, cphi, diurab, era1, refa, refb.
    ///     eo: equation of the origins (i.e., ERA-GST).</returns>
    public static Tuple<iauASTROM, double> Apci13(double date1, double date2) {
        var astrom = new iauASTROM();
        double eo = 0;
        S_Apci13(date1, date2, ref astrom, ref eo);
        return Tuple.Create(astrom, eo); }
    [DllImport(DllFilename, EntryPoint = "iauApci13", CallingConvention = CallingConvention.Cdecl)]
    static extern void S_Apci13(double date1, double date2, ref iauASTROM astrom, ref double eo);

    // Apco():  NOT IMPLEMENTED: earth ephemeris must be supplied by caller, so prefer Apco13(). 

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
    /// <param name="phi">Latitude of observer (radians, geodetic)</param>
    /// <param name="hm">Height above ellipsoid (meters, geodetic)</param>
    /// <param name="xp">Polar motion x coordinate (radians).
    /// Obtained from IERS bulletins, or may be set to zero with little consequence.</param>
    /// <param name="yp">Polar motion y coordinate (radians).
    /// Obtained from IERS bulletins, or may be set to zero with little consequence.</param>
    /// <param name="phpa">Pressure at the observer (hPa = mb). May be estimated from hm as
    /// 1013.25 * exp(-hm / (29.3 * tsl)) where tsl is approx. sea-level temperature Kelvin.</param>
    /// <param name="tc">Ambient temperature at the observer (deg C)</param>
    /// <param name="rh">Relative humidity at the observer (range 0-1)</param>
    /// <param name="w1">Wavelength of observation (micrometers)</param>
    /// <returns> 3-Tuple:
    ///     status: [int] +1 -> dubious year, 0 -> OK, -1 -> unacceptable date.
    ///     astrom: Star−independent astrometry parameters [iauASTROM struct].
    ///             Unchanged struct elements: along, xp1, yp1, sphi, cphi, diurab, era1, refa, refb.
    ///     eo: [double] Equation of the origins (i.e., ERA-GST)</returns>
    public static Tuple<int, iauASTROM, double> Apco13(double utc1, double utc2, double dut1, double elong,
        double phi, double hm, double xp, double yp, double phpa, double tc, double rh, double wl) {
        var astrom = new iauASTROM();
        double eo = 0;
        var status = S_Apco13(utc1, utc2, dut1, elong, phi, hm, xp, yp, phpa, tc, rh, wl, 
            ref astrom, ref eo);
        return Tuple.Create(status, astrom, eo); }
    [DllImport(DllFilename, EntryPoint = "iauApco13", CallingConvention = CallingConvention.Cdecl)]
    static extern int S_Apco13(double utc1, double utc2, double dut1, double elong, 
        double phi, double hm, double xp, double yp, double phpa, double tc, double rh, double wl, 
        ref iauASTROM astrom, ref double eo);

    // Apcs():  NOT IMPLEMENTED: earth ephemeris must be supplied by caller, so prefer Apcs13(). 

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
    /// (Barycentric) axes, and NOT geocentric axes (meters, meters/second)</param>
    /// <returns>Star−independent astrometry parameters, as a iauASTROM struct.
    /// Unchanged struct elements: along, xp1, yp1, sphi, cphi, diurab, era1, refa, refb.</returns>
    public static iauASTROM Apcs13(double date1, double date2, double[,] pv) {
        var astrom = new iauASTROM();
        S_Apcs13(date1, date2, pv, ref astrom);
        return astrom; }
    [DllImport(DllFilename, EntryPoint = "iauApcs13", CallingConvention = CallingConvention.Cdecl)]
    static extern void S_Apcs13(double date1, double date2, double[,] pv, ref iauASTROM astrom);

    // Aper():  NOT IMPLEMENTED: earth ephemeris must be supplied by caller, so prefer Aper13(). 
    
    /// <summary> Sofa.Aper13(): For a TERRESTRIAL OBSERVER, use UT1 and longitude
    /// to update "local" Earth rotation angle only.</summary>
    /// <param name="ut11">UT1 (not UTC) as Julian date, part 1</param>
    /// <param name="ut12">UT1 (not UTC) as Julian date, part 2</param>
    /// <param name="astrom">Astrometry struct iauASTROM, but only the along element is used</param>
    /// <returns>astrom: Star−independent astrometry parameter era1, as the only updated element of
    /// this iauASTROM struct.</returns>
    public static iauASTROM Aper13(double ut11, double ut12, ref iauASTROM astrom) {
        S_Aper13(ut11, ut12, ref astrom);
        return astrom; }
    [DllImport(DllFilename, EntryPoint = "iauAper13", CallingConvention = CallingConvention.Cdecl)]
    static extern void S_Aper13(double ut11, double ut12, ref iauASTROM astrom);

    // Apio():  NOT IMPLEMENTED: earth ephemeris must be supplied by caller, so prefer Apio13(). 

    /// <summary> Sofa.Apio13(): For a TERRESTRIAL OBSERVER, prepare star−independent astrometry parameters
    /// for transformations between CIRS and observed coordinates. The caller supplies UTC, site coordinates,
    /// ambient air conditions and observing wavelength, and SOFA models are used to predict
    /// the Earth ephemeris. The parameters produced by this function are required in the parallax,
    /// light deflection and aberration parts of the astrometric transformation chain.</summary>
    /// <param name="utc1">UTC as quasi-Julian date, part 1. Use Sofa.Dtf2d() to convert from calendar
    /// date and time of day into 2-part quasi-Julian date, as it implements the proper leap-second
    /// ambiguity convention.</param>
    /// <param name="utc2">UTC as quasi-Julian date, part 2</param>
    /// <param name="dut1">UT1-UTC (seconds). Tabulated in IERS bulletins. </param>
    /// <param name="elong">Longitude of observer (radians, WGS84). CAUTION: east positive).</param>
    /// <param name="phi">Latitude of observer (radians, geodetic)</param>
    /// <param name="hm">Height above ellipsoid (meters, geodetic)</param>
    /// <param name="xp">Polar motion x coordinate (radians).
    /// Obtained from IERS bulletins, or may be set to zero with little consequence.</param>
    /// <param name="yp">Polar motion y coordinate (radians).
    /// Obtained from IERS bulletins, or may be set to zero with little consequence.</param>
    /// <param name="phpa">Pressure at the observer (hPa = mb). May be estimated from hm as
    /// 1013.25 * exp(-hm / (29.3 * tsl)) where tsl is approx. sea-level temperature Kelvin.</param>
    /// <param name="tc">Ambient temperature at the observer (deg C)</param>
    /// <param name="rh">Relative humidity at the observer (range 0-1)</param>
    /// <param name="wl">Wavelength of observation (micrometers)</param>
    /// <returns> 2-Tuple:
    ///         status: [int] +1 -> dubious year; 0 -> OK; -1 -> unacceptable date.
    ///         astrom: Star−independent astrometry parameters, as a iauASTROM struct.
    ///                 Struct elements unchanged: pmt, eb, eh, em, v, bm1, bpn.</returns>
    public static Tuple<int, iauASTROM> Apio13(double utc1, double utc2, double dut1, double elong,
        double phi, double hm, double xp, double yp, double phpa, double tc, double rh, double wl) {
        var astrom = new iauASTROM();
        var status = S_Apio13(utc1, utc2, dut1, elong, phi, hm, xp, yp, phpa, tc, rh, wl, ref astrom);
        return Tuple.Create(status, astrom); }
    [DllImport(DllFilename, EntryPoint = "iauApio13", CallingConvention = CallingConvention.Cdecl)]
    static extern int S_Apio13(double utc1, double utc2, double dut1, double elong, double phi, double hm,
        double xp, double yp, double phpa, double tc, double rh, double wl, ref iauASTROM astrom);

    /// <summary>Sofa.Atcc13(): Transform a star’s ICRS catalog entry (epoch J2000.0) into ICRS
    /// astrometric place (update J2000.0 RA,Dec to datetime other than 2000.0).
    /// [NB: this fails for ATLAS refcat2, whose position epoch is NOT 2000.0.)</summary>
    /// <param name="rc">ICRS right ascension at J2000.0 (radians)</param>
    /// <param name="dc">ICRS declination at J2000.0 (radians)</param>
    /// <param name="pr">RA proper motion (radians/year; as dRA/dt, not cos(Dec)*dRA/dt)</param>
    /// <param name="pd">Dec proper motion (radians/year)</param>
    /// <param name="px">Parallax (arcseconds)</param>
    /// <param name="rv">Radial velocity (km/s, +ve if receding)</param>
    /// <param name="date1">Part 1 of 2-part Julian date</param>
    /// <param name="date2">Part 2 of 2-part Julian date</param>
    /// <returns>2-tuple of doubles: RA, Dec (ICRS astrometric, radians).</returns>
    public static Tuple<double, double> Atcc13(double rc, double dc, double pr, double pd,
        double px, double rv, double date1, double date2) {
        double ra = 0, da = 0;
        S_Atcc13(rc, dc, pr, pd, px, rv, date1, date2, ref ra, ref da);
        return Tuple.Create(ra, da); }
    [DllImport(DllFilename, EntryPoint = "iauAtcc13", CallingConvention = CallingConvention.Cdecl)]
    static extern void S_Atcc13(double rc, double dc, double pr, double pd, double px, double rv,
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
    /// by another SOFA function.</param>
    /// <returns>2-tuple:  ra: ICRS astrometric Right Ascension (radians).
    ///                    da: ICRS astrometric Declination (radians).</returns>
    public static Tuple<double, double> Atccq(double rc, double dc, double pr, double pd,
        double px, double rv, iauASTROM astrom) {
        double ra = 0, da = 0;
        S_Atccq(rc, dc, pr, pd, px, rv, astrom, ref ra, ref da);
        return Tuple.Create(ra, da); }
    [DllImport(DllFilename, EntryPoint = "iauAtccq", CallingConvention = CallingConvention.Cdecl)]
    static extern void S_Atccq(double rc, double dc, double pr, double pd, double px, double rv,
        iauASTROM astrom, ref double ra, ref double da);

    /// <summary>Sofa.Atci13(): Transform ICRS star data, epoch J2000.0, to CIRS.
    /// NB: SOFA param namings are reversed: rc & dc are ICRS, not CIRS; ri & di are CIRS.</summary>
    /// <param name="rc">ICRS right ascension at J2000.0 (radians)</param>
    /// <param name="dc">ICRS declination at J2000.0 (radians)</param>
    /// <param name="pr">RA proper motion (radians/year; as dRA/dt, not cos(Dec)*dRA/dt)</param>
    /// <param name="pd">Dec proper motion (radians/year)</param>
    /// <param name="px">Parallax (arcseconds)</param>
    /// <param name="rv">Radial velocity (km/s, +ve if receding)</param>
    /// <param name="date1">Part 1 of 2-part Julian date</param>
    /// <param name="date2">Part 2 of 2-part Julian date</param>
    /// <returns>3-tuple of doubles: RA (ICRS astrometric),
    ///                              Dec (ICRS astrometric),
    ///                              eo (equation of origins, ERA-GST).</returns>
    public static Tuple<double, double, double> Atci13(double rc, double dc, double pr, double pd,
        double px, double rv, double date1, double date2) {
        double ri = 0, di = 0, eo = 0;
        S_Atci13(rc, dc, pr, pd, px, rv, date1, date2, ref ri, ref di, ref eo);
        return Tuple.Create(ri, di, eo); }
    [DllImport(DllFilename, EntryPoint = "iauAtci13", CallingConvention = CallingConvention.Cdecl)]
    static extern void S_Atci13(double rc, double dc, double pr, double pd, double px, double rv,
        double date1, double date2, ref double ri, ref double di, ref double eo);
    
    /// <summary> Sofa.Atciq(): Quick ICRS (epoch J2000.0) to CIRS transformation, given precomputed
    /// star−independent astrometry parameters.
    ///
    /// Use of this function is appropriate when efficiency is important and where many star positions
    /// are to be transformed for one date. The star−independent parameters can be obtained by
    /// calling one of the functions iauApci[13], iauApcg[13], iauApco[13] or iauApcs[13].
    ///
    /// If the parallax and proper motions are zero the iauAtciqz function can be used instead.</summary>
    /// <param name="rc">ICRS Right Ascension (radians, J2000.0)</param>
    /// <param name="dc">ICRS Declination (radians, J2000.0)</param>
    /// <param name="pr">Proper motion in Right Ascension (radians/year)</param>
    /// <param name="pd">Proper motion in Declination (radians/year)</param>
    /// <param name="px">Parallax (arcseconds)</param>
    /// <param name="rv">Radial velocity (km/s, receding is positive)</param>
    /// <param name="astrom">Star-independent astrometry parameters, as a iauASTROM struct generated
    /// by another SOFA function</param>
    /// <returns>2-tuple:  ri: CIRS astrometric Right Ascension (radians)
    ///                    di: CIRS astrometric Declination (radians)</returns>
    public static Tuple<double, double> Atciq(double rc, double dc, double pr, double pd,
        double px, double rv, iauASTROM astrom) {
        double ri = 0, di = 0;
        S_Atciq(rc, dc, pr, pd, px, rv, astrom, ref ri, ref di);
        return Tuple.Create(ri, di); }
    [DllImport(DllFilename, EntryPoint = "iauAtciq", CallingConvention = CallingConvention.Cdecl)]
    static extern void S_Atciq(double rc, double dc, double pr, double pd, double px, double rv,
        iauASTROM astrom, ref double ri, ref double di);

    // Atciqn():  NOT IMPLEMENTED: light-deflecting bodies are of low relevance to amateur astronomy. 

    /// <summary> Sofa.Atciqz(): Quick ICRS (epoch J2000.0) to CIRS transformation, given precomputed
    /// star−independent astrometry parameters, and assuming zero parallax and proper motion.
    ///
    /// Use of this function is appropriate when efficiency is important and where many star positions
    /// are to be transformed for one date. The star−independent parameters can be obtained by
    /// calling one of the functions iauApci[13], iauApcg[13], iauApco[13] or iauApcs[13].
    ///
    /// The corresponding function for the case of non−zero parallax and proper motion
    /// is iauAtciq.</summary>
    /// <param name="rc">ICRS Right Ascension (radians, J2000.0)</param>
    /// <param name="dc">ICRS Declination (radians, J2000.0)</param>
    /// <param name="astrom">Star-independent astrometry parameters, as a iauASTROM struct generated
    /// by another SOFA function</param>
    /// <returns>2-tuple:  ri: CIRS astrometric Right Ascension (radians)
    ///                    di: CIRS astrometric Declination (radians)</returns>
    public static Tuple<double, double> Atciqz(double rc, double dc, iauASTROM astrom) {
        double ri = 0, di = 0;
        S_Atciqz(rc, dc, astrom, ref ri, ref di);
        return Tuple.Create(ri, di); }
    [DllImport(DllFilename, EntryPoint = "iauAtciqz", CallingConvention = CallingConvention.Cdecl)]
    static extern void S_Atciqz(double rc, double dc, iauASTROM astrom, ref double ri, ref double di);

    /// <summary>Sofa.Atco13(): ICRS RA,Dec to observed place. The caller supplies UTC,
    /// site coordinates, ambient air conditions and observing wavelength.</summary>
    /// <param name="rc">ICRS right ascension at J2000.0 (radians)</param>
    /// <param name="dc">ICRS declination at J2000.0 (radians)</param>
    /// <param name="pr">RA proper motion (radians/year; as dRA/dt, not cos(Dec)*dRA/dt)</param>
    /// <param name="pd">Dec proper motion (radians/year)</param>
    /// <param name="px">Parallax (arcseconds)</param>
    /// <param name="rv">Radial velocity (km/s, +ve if receding)</param>
    /// <param name="utc1">Part 1 of 2-part UTC quasi-Julian date</param>
    /// <param name="utc2">Part 2 of 2-part UTC quasi-Julian date</param>
    /// <param name="dut1">UT1−UTC (seconds, Note 5)</param>
    /// <param name="elong">Longitude (radians, WGS84, east +ve)</param>
    /// <param name="phi">Latitude (radians, geodetic, WGS84)</param>
    /// <param name="hm">Height above ellipsoid (meters, geodetic, WGS84)</param>
    /// <param name="xp">Polar motion coordinate (radians, from IERS bulletins, often ~zero)</param>
    /// <param name="yp">Polar motion coordinate (radians, from IERS bulletins, often ~zero)</param>
    /// <param name="phpa">Pressure at the observer (hPa = mB)</param>
    /// <param name="tc">Ambient temperature at the observer (degrees C)</param>
    /// <param name="rh">Relative humidity at the observer (in range [0-1])</param>
    /// <param name="wl">Observation wavelength (micrometers)</param> 
    /// <returns>7-tuple: status [int]: +1 -> dubious year; 0 -> OK; -1 -> unacceptable date.
    ///                   aob: observed azimuth (radians),
    ///                   zob: observed zenith distance NB: NOT altitude (radians),
    ///                   hob: observed hour angle (radians),
    ///                   dob: observed declination (radians),
    ///                   rob: observed right ascension (radians, CIO-based),
    ///                   eo:  equation of the origins (ERA, GST).</returns>
    public static Tuple<int, double, double, double, double, double, double> Atco13(double rc, double dc,
        double pr, double pd, double px, double rv, double utc1, double utc2, double dut1,
        double elong, double phi, double hm, double xp, double yp, double phpa,
        double tc, double rh, double wl) {
        double aob = 0, zob = 0, hob = 0, dob = 0, rob = 0, eo = 0;
        var status = S_Atco13(rc, dc, pr, pd, px, rv, utc1, utc2, dut1, elong, phi, hm, xp, yp, phpa, 
            tc, rh, wl, ref aob, ref zob, ref hob, ref dob, ref rob, ref eo);
        return Tuple.Create(status, aob, zob, hob, dob, rob, eo); }
    [DllImport(DllFilename, EntryPoint = "iauAtco13", CallingConvention = CallingConvention.Cdecl)]
    static extern int S_Atco13(double rc, double dc, double pr, double pd, double px, double rv,
        double utc1, double utc2, double dut1, double elong, double phi, double hm, double xp, double yp,
        double phpa, double tc, double rh, double wl,
        ref double aob, ref double zob, ref double hob, ref double dob, ref double rob, ref double eo);
    
    /// <summary>Sofa.Atic13(): Transform star RA,Dec from geocentric CIRS to ICRS astrometric.</summary>
    /// <param name="ri">CIRS geocentric Right Ascension (radians)</param>
    /// <param name="di">CIRS geocentric Declination (radians)</param>
    /// <param name="date1">Part 1 of 2-part Julian date</param>
    /// <param name="date2">Part 2 of 2-part Julian date</param>
    /// <returns>3-tuple of doubles: rc: Right Ascension (ICRS astrometric),
    ///                              dc: Declination (ICRS astrometric),
    ///                              eo: Equation of origins, ERA-GST.</returns>
    public static Tuple<double, double, double> Atic13(double ri, double di, double date1, double date2) {
        double rc = 0, dc = 0, eo = 0;
        S_Atic13(ri, di, date1, date2, ref rc, ref dc, ref eo);
        return Tuple.Create(rc, dc, eo); }
    [DllImport(DllFilename, EntryPoint = "iauAtic13", CallingConvention = CallingConvention.Cdecl)]
    static extern void S_Atic13(double ri, double di, double date1, double date2,
        ref double rc, ref double dc, ref double eo);

    /// <summary> Sofa.Aticq(): Quick transform from CIRS RA,Dec to ICRS astrometric place,
    /// given star-independent astrometry parameters.</summary>
    /// <param name="ri">CIRS geocentric RA (radians)</param>
    /// <param name="di">CIRS geocentric Dec (radians)</param>
    /// <param name="astrom">Star-independent astrometry parameters, as a iauASTROM struct generated
    /// by another SOFA function</param>
    /// <returns>2-tuple of doubles: rc, dc: RA, Dec (ICRS astrometric).</returns>
    public static Tuple<double, double> Aticq(double ri, double di, iauASTROM astrom) {
        double rc = 0, dc = 0;
        S_Aticq(ri, di, astrom, ref rc, ref dc);
        return Tuple.Create(rc, dc); }
    [DllImport(DllFilename, EntryPoint = "iauAticq", CallingConvention = CallingConvention.Cdecl)]
    static extern void S_Aticq(double ri, double di, iauASTROM astrom, ref double rc, ref double dc);

    // Aticqn():  NOT IMPLEMENTED: light-deflecting bodies are of low relevance to amateur astronomy. 
    
    /// <summary> Sofa.Atio13(): Transform CIRS RA,Dec to observed place. The caller supplies UTC,
    /// site coordinates, ambient air conditions and observing wavelength.</summary>
    /// <param name="ri">CIRS geocentric RA (radians)</param>
    /// <param name="di">CIRS geocentric Dec (radians)</param>
    /// <param name="utc1">Part 1 of 2-part UTC quasi-Julian date</param>
    /// <param name="utc2">Part 2 of 2-part UTC quasi-Julian date</param>
    /// <param name="dut1">UT1−UTC (seconds, Note 5)</param>
    /// <param name="elong">Longitude (radians, WGS84, east +ve)</param>
    /// <param name="phi">Latitude (radians, geodetic, WGS84)</param>
    /// <param name="hm">Height above ellipsoid (meters, geodetic, WGS84)</param>
    /// <param name="xp">Polar motion coordinate (radians, from IERS bulletins, often ~zero)</param>
    /// <param name="yp">Polar motion coordinate (radians, from IERS bulletins, often ~zero)</param>
    /// <param name="phpa">Pressure at the observer (hPa = mB)</param>
    /// <param name="tc">Ambient temperature at the observer (degrees C)</param>
    /// <param name="rh">Relative humidity at the observer (in range [0-1])</param>
    /// <param name="wl">Observation wavelength (micrometers)</param> 
    /// <returns>6-tuple of doubles: status [int]: +1 -> dubious year; 0 -> OK; -1 -> unacceptable date,
    ///                              aob: observed azimuth (radians),
    ///                              zob: observed zenith distance (radians),
    ///                              hob: observed hour angle (radians),
    ///                              dob: observed declination (radians),
    ///                              rob: observed right ascension (radians, CIO-based).</returns>
    public static Tuple<int, double, double, double, double, double> Atio13(double ri, double di,
        double utc1, double utc2, double dut1, double elong, double phi, double hm, double xp, double yp,
        double phpa, double tc, double rh, double wl) {
        double aob = 0, zob = 0, hob = 0, dob = 0, rob = 0;
        var status = S_Atio13(ri, di, utc1, utc2, dut1, elong, phi, hm, xp, yp, phpa, tc, rh, wl,
            ref aob, ref zob, ref hob, ref dob, ref rob);
        return Tuple.Create(status, aob, zob, hob, dob, rob); }
    [DllImport(DllFilename, EntryPoint = "iauAtio13", CallingConvention = CallingConvention.Cdecl)]
    static extern int S_Atio13(double ri, double di, double utc1, double utc2,
        double dut1, double elong, double phi, double hm, double xp, double yp,
        double phpa, double tc, double rh, double wl,
        ref double aob, ref double zob, ref double hob, ref double dob, ref double rob);

    /// <summary> Sofa.Atioq(): Quick transform from CIRS to Observed place.
    /// Use this function when efficiency is important and where many star positions are all
    /// to be transformed for one date. The star−independent astrometry parameters can be obtained by
    /// calling iauApio[13] or iauApco[13].</summary>
    /// <param name="ri">CIRS Right Ascension (radians)</param>
    /// <param name="di">CIRS Declination (radians)</param>
    /// <param name="astrom">Star-independent astrometry parameters, as a iauASTROM struct generated
    /// by another SOFA function</param>
    /// <returns>5-tuple of doubles: aob: observed azimuth (radians, N=0, E=pi/2),
    ///                              zob: observed zenith distance (radians, zenith=0),
    ///                              hob: observed hour angle (radians),
    ///                              dob: observed declination (radians),
    ///                              rob: observed right ascension (radians, CIO-based)</returns>
    public static Tuple<double, double, double, double, double> Atioq(double ri, double di,
        iauASTROM astrom) {
        double aob = 0, zob = 0, hob = 0, dob = 0, rob = 0;
        S_Atioq(ri, di, astrom, ref aob, ref zob, ref hob, ref dob, ref rob);
        return Tuple.Create(aob, zob, hob, dob, rob); }
    [DllImport(DllFilename, EntryPoint = "iauAtioq", CallingConvention = CallingConvention.Cdecl)]
    static extern void S_Atioq(double ri, double di, iauASTROM astrom,
        ref double aob, ref double zob, ref double hob, ref double dob, ref double rob);

    /// <summary>Sofa.Atoc13(): Transform observed place at a ground-based site to ICRS
    /// astrometric RA,Dec. The caller supplies UTC, site coordinates,
    /// ambient air conditions, and observing wavelength.</summary>
    /// <param name="type">type of coordinates (string of length one, case-insensitive):
    ///                    "R" -> ob1=RA, ob2=Dec,
    ///                    "H" -> ob1=Hour Angle, ob2=Dec, or
    ///                    "A" -> ob1=Azimuth (N=0, E=90deg), ob2=zenith distance).</param>
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
    /// <returns>3-tuple: status [int]: +1 -> dubious year; 0 -> OK; -1 -> unacceptable date,
    ///                   rc: Right Ascension (ICRS astrometric),
    ///                   dc: Declination (ICRS astrometric).</returns>.
    public static Tuple<int, double, double> Atoc13(string type, double ob1, double ob2, double utc1,
        double utc2,
        double dut1, double elong, double phi, double hm, double xp, double yp,
        double phpa, double tc, double rh, double wl) {
        double rc = 0, dc = 0;
        var status = S_Atoc13(type, ob1, ob2, utc1, utc2, dut1, elong, phi, hm, xp, yp, phpa, tc, rh, wl,
            ref rc, ref dc);
        return Tuple.Create(status, rc, dc); }
    [DllImport(DllFilename, EntryPoint = "iauAtoc13", CallingConvention = CallingConvention.Cdecl)]
    static extern int S_Atoc13(string type, double ob1, double ob2, double utc1, double utc2,
        double dut1, double elong, double phi, double hm, double xp, double yp, double phpa, double tc,
        double rh, double wl, ref double rc, ref double dc);

    /// <summary>Sofa.Atoi13(): Transform observed place to CIRS. The caller supplies UTC,
    /// site coordinates, ambient air conditions and observing wavelength.</summary>
    /// <param name="type">type of coordinates (string of length one, case-insensitive):
    ///                    "R" -> ob1=RA, ob2=Dec,
    ///                    "H" -> ob1=Hour Angle, ob2=Dec, or
    ///                    "A" -> ob1=Azimuth (N=0, E=90deg), ob2=zenith distance).</param>
    /// <param name="ob1">observed RA, Hour Angle, or Azimuth (N=0, E+) (radians)</param>
    /// <param name="ob2">observed Zenith Distance or Declination (radians)</param>
    /// <param name="utc1">Part 1 of 2-part UTC quasi-Julian date</param>
    /// <param name="utc2">Part 2 of 2-part UTC quasi-Julian date</param>
    /// <param name="dut1">UT1−UTC (seconds, Note 5)</param>
    /// <param name="elong">Longitude (radians, WGS84, east +ve)</param>
    /// <param name="phi">Latitude (radians, geodetic, WGS84)</param>
    /// <param name="hm">Height above ellipsoid (meters, geodetic, WGS84)</param>
    /// <param name="xp">Polar motion coordinate (radians, from IERS bulletins, often ~zero)</param>
    /// <param name="yp">Polar motion coordinate (radians, from IERS bulletins, often ~zero)</param>
    /// <param name="phpa">Pressure at the observer (hPa = mB)</param>
    /// <param name="tc">Ambient temperature at the observer (degrees C)</param>
    /// <param name="rh">Relative humidity at the observer (in range [0-1])</param>
    /// <param name="wl">Observation wavelength (micrometers)</param> 
    /// <returns>3-tuple: status [int]: +1 -> dubious year; 0 -> OK; -1 -> unacceptable date,
    ///                   ri: Right Ascension (radians, CIRS, CIO-based),
    ///                   di: Declination (radians, CIRS, CIO-based).</returns>.
    public static Tuple<int, double, double> Atoi13(string type, double ob1, double ob2, 
        double utc1, double utc2, double dut1, double elong, double phi, double hm, double xp, double yp,
        double phpa, double tc, double rh, double wl) {
        double ri = 0, di = 0;
        var status = S_Atoi13(type, ob1, ob2, utc1, utc2, dut1, elong, phi, hm, xp, yp, phpa, 
            tc, rh, wl, ref ri, ref di);
        return Tuple.Create(status, ri, di); }
    [DllImport(DllFilename, EntryPoint = "iauAtoi13", CallingConvention = CallingConvention.Cdecl)]
    static extern int S_Atoi13(string type, double ob1, double ob2, double utc1, double utc2,
        double dut1, double elong, double phi, double hm, double xp, double yp, double phpa, double tc,
        double rh, double wl, ref double ri, ref double di);
    
    /// <summary> Sofa.Atoiq(): Quick transform from Observed place to CIRS,
    /// given the star-independent astrometry parameters.
    /// Use this function when efficiency is important and where many star positions are all to be
    /// transformed for one date. The star−independent astrometry parameters can be obtained by
    /// calling iauApio[13] or iauApco[13].</summary>
    /// <param name="type">type of coordinates (string of length one, case-insensitive):
    ///                    "R" -> ob1=RA, ob2=Dec,
    ///                    "H" -> ob1=Hour Angle, ob2=Dec, or
    ///                    "A" -> ob1=Azimuth (N=0, E=90deg), ob2=zenith distance).</param>
    /// <param name="ob1">observed RA, Hour Angle, or Azimuth (N=0, E+) (radians)</param>
    /// <param name="ob2">observed Zenith Distance or Declination (radians)</param>
    /// <param name="astrom">Star-independent astrometry parameters, as a iauASTROM struct generated
    /// by another SOFA function</param>
    /// <returns>2-tuple of doubles: ri: Right Ascension (radians, CIRS, CIO-based),
    ///                              di: Declination (radians, CIRS, CIO-based).</returns>.
    public static Tuple<double, double> Atoiq(string type, double ob1, double ob2, iauASTROM astrom) {
        double ri = 0, di = 0;
        S_Atoiq(type, ob1, ob2, astrom, ref ri, ref di);
        return Tuple.Create(ri, di); }
    [DllImport(DllFilename, EntryPoint = "iauAtoiq", CallingConvention = CallingConvention.Cdecl)]
    static extern void S_Atoiq(string type, double ob1, double ob2, iauASTROM astrom,
        ref double ri, ref double di);

    // Ld():    NOT IMPLEMENTED: light deflection, which is of low relevance to amateur astronomy.
    // Ldn():   NOT IMPLEMENTED: light deflection, which is of low relevance to amateur astronomy.
    // Lsun():  NOT IMPLEMENTED: light deflection, which is of low relevance to amateur astronomy.
    // Pmpx():  NOT IMPLEMENTED: Proper motion -> BCRS unit vector, of low relevance to amateur astronomy.
    
    /// <summary>Sofa.Pmsafe(): Star proper motion: update star catalog data for space motion, with
    ///                         special handling to handle the zero parallax case.</summary>
    /// <param name="ra1">Right Ascension, before (radians)</param>
    /// <param name="dec1">Declination, before (radians)</param>
    /// <param name="pmr1">RA Proper motion, before (radians/year)</param>
    /// <param name="pmd1">Declination Proper motion, before (radians/year)</param>
    /// <param name="px1">Parallax, before (arcseconds)</param>
    /// <param name="rv1">Radial velocity, before (km/s, +v -> receding)</param>
    /// <param name="ep1a">Part A of 2-part Epoch TDB quasi-JD, before</param>
    /// <param name="ep1b">Part B of 2-part Epoch TDB quasi-JD, before</param>
    /// <param name="ep2a">Part A of 2-part Epoch TDB quasi-JD, after</param>
    /// <param name="ep2b">Part B of 2-part Epoch TDB quasi-JD, after</param>
    /// <returns>7-tuple:  status:  [int] -1 -> system error (should never occur)
    ///                                    0 -> OK
    ///                                    1 -> distance overridden (too small or negative)
    ///                                    2 -> excessive velocity (then arbitrarily set to zero)
    ///                                    4 -> solution did not converge
    ///                                    other -> multiple warnings (binary logical OR of the above)
    ///                    ra2:     Right Ascension, after (radians)
    ///                    dec2:    Declination, after (radians)
    ///                    pmr2:    RA proper motion, after (radians/year)
    ///                    pmd2:    Declination proper motion, after (radians/year)
    ///                    px2:     Parallax, after (arcseconds)
    ///                    rv2:     Radial velocity, after (km/s, +v -> receding). </returns>.
    public static Tuple<int, double, double, double, double, double, double> Pmsafe(
        double ra1, double dec1, double pmr1, double pmd1, double px1, double rv1,
        double ep1a, double ep1b, double ep2a, double ep2b) {
        double ra2 = 0, dec2 = 0, pmr2 = 0, pmd2 = 0, px2 = 0, rv2 = 0;
        var status = S_Pmsafe(ra1, dec1, pmr1, pmd1, px1, rv1, ep1a, ep1b, ep2a, ep2b,
            ref ra2, ref dec2, ref pmr2, ref pmd2, ref px2, ref rv2);
        return Tuple.Create(status, ra2, dec2, pmr2, pmd2, px2, rv2); }
    [DllImport(DllFilename, EntryPoint = "iauPmsafe", CallingConvention = CallingConvention.Cdecl)]
    static extern int S_Pmsafe(double ra1, double dec1, double pmr1, double pmd1,
        double px1, double rv1, double ep1a, double ep1b, double ep2a, double ep2b,
        ref double ra2, ref double dec2, ref double pmr2, ref double pmd2, ref double px2, ref double rv2);

    // Pvstar():   NOT IMPLEMENTED: Assumes velocity of observed body is known without catalog.
    
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
        var pv = new double[2, 3] {{0, 0, 0}, {0, 0, 0}};
        S_Pvtob(elong, phi, hm, xp, yp, sp, theta, pv);
        return pv; }
    [DllImport(DllFilename, EntryPoint = "iauPvtob", CallingConvention = CallingConvention.Cdecl)]
    static extern void S_Pvtob(double elong, double phi, double hm, double xp, double yp,
        double sp, double theta, [In, Out] double[,] pv);

    /// <summary>Sofa.Refco(): Determine the constants A and B in the atmospheric refraction model
    /// dZ = A tan Z + B tan^3 Z.</summary>
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
    static extern void S_Refco(double phpa, double tc, double rh, double wl,
        ref double refa, ref double refb);

    /// <summary>Sofa.Starpm(): Star proper motion: update (for new time) star catalog data
    ///                         for space motion.
    /// NB: not safe for zero parallax and/or radial velocity; when in doubt, use Sofa.Pmsafe().</summary>
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
    /// <returns>7-tuple:  status:  [int] -1 -> system error (should never occur)
    ///                                    0 -> OK
    ///                                    1 -> distance overridden (too small or negative)
    ///                                    2 -> excessive velocity (then arbitrarily set to zero)
    ///                                    4 -> solution did not converge
    ///                                    other -> multiple warnings (binary logical OR of the above)
    ///                    ra2:  Right Ascension, after (radians)
    ///                    dec2: Declination, after (radians)
    ///                    pmr2: RA proper motion, after (radians/year)
    ///                    pmd2: Declination proper motion, after (radians/year)
    ///                    px2:  Parallax, after (arcseconds)
    ///                    rv2:  Radial velocity, after (km/s, positive=receding).</returns>.
    public static Tuple<int, double, double, double, double, double, double> 
        Starpm(double ra1, double dec1, double pmr1, double pmd1, double px1, double rv1,
        double ep1a, double ep1b, double ep2a, double ep2b) {
        if (px1 == 0.0 && rv1 == 0.0) {
            throw new ArgumentException(
                "Use Sofa.Pmsafe() for zero parallax and/or zero radial velocity."); }
        double ra2 = 0, dec2 = 0, pmr2 = 0, pmd2 = 0, px2 = 0, rv2 = 0;
        var status = S_Starpm(ra1, dec1, pmr1, pmd1, px1, rv1, ep1a, ep1b, ep2a, ep2b,
            ref ra2, ref dec2, ref pmr2, ref pmd2, ref px2, ref rv2);
        return Tuple.Create(status, ra2, dec2, pmr2, pmd2, px2, rv2); }
    [DllImport(DllFilename, EntryPoint = "iauStarpm", CallingConvention = CallingConvention.Cdecl)]
    static extern int S_Starpm(double ra1, double dec1, double pmr1, double pmd1, double px1, double rv1,
        double ep1a, double ep1b, double ep2a, double ep2b,
        ref double ra2, ref double dec2, ref double pmr2, ref double pmd2, ref double px2, ref double rv2);

    // Starpv():   NOT IMPLEMENTED: Produces position/velocity vector,
    //             usually of low relevance to amateur astronomy (may implement later).

    #endregion  //  Astrometry


    #region "Astronomy: Calendars"
    /*######### Astronomy: Calendars #####################################################################*/
    
    /// <summary>Sofa.Cal2jd(): Gregorian Calendar to Julian Date.</summary>
    /// <param name="iy">Year in Gregorian calendar</param>
    /// <param name="im">Month in Gregorian calendar</param>
    /// <param name="id">Day in Gregorian calendar</param>
    /// <returns>3-tuple: status [int]:  0 -> OK,
    ///                                 -1 -> bad year (JD not computed),
    ///                                 -2 -> bad month (JD not computed),
    ///                                 -3 -> bad day (JD computed but probably wrong).
    ///                   djm0: MJD zero−point, always 2400000.5,
    ///                   djm:  Modified Julian Date for 0 hrs. </returns>
    public static Tuple<int, double, double> Cal2jd(int iy, int im, int id) {
        double djm0 = 0, djm = 0;
        var status = S_Cal2jd(iy, im, id, ref djm0, ref djm);
        return Tuple.Create(status, djm0, djm); }
    [DllImport(DllFilename, EntryPoint = "iauCal2jd", CallingConvention = CallingConvention.Cdecl)]
    static extern int S_Cal2jd(int iy, int im, int id, ref double djm0, ref double djm);

    
    /// <summary>Sofa.Jd2cal(): Julian Date to Gregorian year, month, day, and fraction of a day.</summary>
    /// <param name="dj1">Part 1 of 2-part Julian Date</param>
    /// <param name="dj2">Part 2 of 2-part Julian Date</param>
    /// <returns>5-tuple: status [int]: 0 -> OK, -1 -> unacceptable date (before -4900 March 1).
    ///                   iy, im, id [each int]: year, month, and day,
    ///                   fd [double]: fraction of day.</returns>
    public static Tuple<int, int, int, int, double> Jd2cal(double dj1, double dj2) {
        int iy = 0, im = 0, id = 0;
        double fd = 0;
        var status = S_Jd2cal(dj1, dj2, ref iy, ref im, ref id, ref fd);
        return Tuple.Create(status, iy, im, id, fd); }
    [DllImport(DllFilename, EntryPoint = "iauJd2cal", CallingConvention = CallingConvention.Cdecl)]
    static extern int S_Jd2cal(double dj1, double dj2, ref int iy, ref int im, ref int id, ref double fd);


    /// <summary>Sofa.Jdcalf(): Julian Date to Gregorian Calendar, expressed in a form convenient
    /// for formatting messages: rounded to a specified precision.</summary>
    /// <param name="ndp">Number of decimal places for days in fraction returned
    /// (9 or less for 32-bit integers.).</param>
    /// <param name="dj1">Part 1 of 2-part Julian Date</param>
    /// <param name="dj2">Part 2 of 2-part Julian Date</param>
    /// <returns>iymdf: int[4], where integers are year, month, day, and fractional representation to
    /// specified precision.</returns>
    public static Tuple<int, int[]> Jdcalf(int ndp, double dj1, double dj2) {
        var iymdf = new int[4] {0, 0, 0, 0};
        var status = S_Jdcalf(ndp, dj1, dj2, iymdf);
        return Tuple.Create(status, iymdf); }
    [DllImport(DllFilename, EntryPoint = "iauJdcalf", CallingConvention = CallingConvention.Cdecl)]
    static extern int S_Jdcalf(int ndp, double dj1, double dj2, [In, Out] int[] iymdf);

    // Epb():      NOT IMPLEMENTED: Besselian Epoch of low relevance to amateur astronomy.
    // Epj():      NOT IMPLEMENTED: Julian Epoch of low relevance to amateur astronomy.
    // Epb2jd():   NOT IMPLEMENTED: Besselian Epoch of low relevance to amateur astronomy.
    // Epj2jd():   NOT IMPLEMENTED: Julian Epoch of low relevance to amateur astronomy.
    
    #endregion //  Calendars


    #region "Astronomy: Time Scales"
    /*######### Astronomy: Time Scales ###################################################################*/

    // From SOFA documentation:
    //-----------------------------------------------------------------
    //   TAI                    <− physically realized
    //    |
    //   offset                 <− observed (nominally +32.184s)
    //    |
    //   TT                     <− terrestrial time
    //    |
    //   rate adjustment (L_G)  <− definition of TT
    //    |
    //   TCG                    <− time scale for GCRS
    //    |
    //   "periodic" terms       <− iauDtdb is an implementation
    //    |
    //   rate adjustment (L_C)  <− function of solar−system ephemeris
    //    |
    //   TCB                    <− time scale for BCRS
    //    |
    //   rate adjustment (−L_B) <− definition of TDB
    //    |
    //   TDB                    <− TCB scaled to track TT
    //    |
    //   "periodic" terms       <− −iauDtdb is an approximation
    //    |
    //   TT                     <− terrestrial time
    //-----------------------------------------------------------------
    // Adopted values for the various constants can be found in the
    //     IERS Conventions (McCarthy & Petit 2003).
    //-----------------------------------------------------------------


    /// <summary>Sofa.D2dtf(): Format for output a 2-part Julian Date
    /// (or in the case of UTC a quasi-JD form that includes special provision for leap seconds).</summary>
    /// <param name="scale">Time scale ID: only "UTC" is significant, enabling handling of leap seconds.
    /// </param>
    /// <param name="ndp">resolution</param>
    /// <param name="d1">first part of time as a 2-part Julian Date</param>
    /// <param name="d2">second part of time as a 2-part Julian Date</param>
    /// <returns>4-tuple: (int year, int month, int day,
    ///                    int array [hours, minutes, seconds, fraction])</returns>
    public static Tuple<int, int, int, int, int[]> D2dtf(string scale, int ndp, double d1, double d2) {
        int iy = 0, im = 0, id = 0;
        var ihmsf = new Int32[4] {0, 0, 0, 0};
        var status = S_D2dtf(scale, ndp, d1, d2, ref iy, ref im, ref id, ihmsf);
        return Tuple.Create(status, iy, im, id, ihmsf); }
    [DllImport(DllFilename, EntryPoint = "iauD2dtf", CallingConvention = CallingConvention.Cdecl)]
    static extern int S_D2dtf(string scale, int ndp, double d1, double d2,
        ref int iy, ref int im, ref int id, [In, Out] int[] ihmsf);


    /// <summary>Sofa.Dat(): For a given UTC date, calculate Delta(AT) = TAI−UTC</summary>
    /// <param name="iy">Year (UTC)"</param>
    /// <param name="im">Month (UTC)</param>
    /// <param name="id">Day (UTC)</param>
    /// <param name="fd">Fraction of day (UTC)</param>
    /// <returns>2-Tuple:  status [int]:  1 -> dubious year (before 1960-01-01)
    ///                                   0 -> OK
    ///                                  -1 -> bad year
    ///                                  -2 -> bad month
    ///                                  -3 -> bad day
    ///                                  -4 -> bad fraction of day, not in [0,1]
    ///                                  -5 -> internal error (rare)
    ///                    deltat: TAI minus UTC (seconds)</returns>
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
    public static Tuple<int, double> Dat(int iy, int im, int id, double fd) {
        double deltat = 0;
        var status = S_Dat(iy, im, id, fd, ref deltat);
        return Tuple.Create(status, deltat); }
    [DllImport(DllFilename, EntryPoint = "iauDat", CallingConvention = CallingConvention.Cdecl)]
    static extern int S_Dat(int iy, int im, int id, double fd, ref double deltat);


    /// <summary>Sofa.Dtdb(): An approximation to TDB−TT, the difference between barycentric
    /// dynamical time and terrestrial time, for an observer on the Earth.</summary>
    /// <param name="date1">Part 1 of 2-part date, TDB</param>
    /// <param name="date2">Part 2 of 2-part date, TDB</param>
    /// <param name="ut">Universal time (UT1, fraction of one day)</param>
    /// <param name="elong">Longitude (radians, east positive)</param>
    /// <param name="u">Distance from Earth spin axis (km)</param>
    /// <param name="v">Distance north of equatorial plane (km) </param>
    /// <returns>double: TDB minus TT (seconds)</returns>
    public static double Dtdb(double date1, double date2, double ut, double elong, double u, double v) {
        var diff = S_Dtdb(date1, date2, ut, elong, u, v);
        return diff; }
    [DllImport(DllFilename, EntryPoint = "iauDtdb", CallingConvention = CallingConvention.Cdecl)]
    static extern double S_Dtdb(double date1, double date2, double ut, double elong, double u, double v);


    /// <summary>Sofa.Dtf2d(): Encode date and time fields into 2−part Julian Date (or in the case
    /// of UTC a quasi−JD form that includes special provision for leap seconds).</summary>
    /// <param name="scale">Timescale ID, typically "UTC"</param>
    /// <param name="iy">Year in Gregorian calendar</param>
    /// <param name="im">Month in Gregorian calendar</param>
    /// <param name="id">Day in Gregorian calendar</param>
    /// <param name="ihr">Hour, UTC</param>
    /// <param name="imn">Minute, UTC</param>
    /// <param name="sec">Seconds, UTC</param>
    /// <returns>3-Tuple:
    ///     status [int]:          +3 -> both +2 and +1 error
    ///                            +2 -> time is after end of day
    ///                            +1 -> dubious year
    ///                             0 -> OK
    ///        -1, -2, -3, -4, -5, -6 -> bad year, month, day, hour, minute, second
    ///     d1, d2, a 2−part Julian Date</returns>
    public static Tuple<int, double, double> Dtf2d(string scale, int iy, int im, int id,
        int ihr, int imn, double sec) {
        double d1 = 0, d2 = 0;
        var status = S_Dtf2d(scale, iy, im, id, ihr, imn, sec, ref d1, ref d2);
        return Tuple.Create(status, d1, d2); }
    [DllImport(DllFilename, EntryPoint = "iauDtf2d", CallingConvention = CallingConvention.Cdecl)]
    static extern int S_Dtf2d(string scale, int iy, int im, int id, int ihr, int imn, double sec,
        ref double d1, ref double d2);


    /// <summary>Sofa.Taitt(): Time scale transformation: International Atomic Time, TAI,
    /// to Terrestrial Time, TT.</summary>
    /// <param name="tai1">Part 1 of 2-part Julian Date, TAI</param>
    /// <param name="tai2">Part 2 of 2-part Julian Date, TAI</param>
    /// <returns>3-Tuple:  status [int]:  0 -> OK
    ///                    tt1: Part 1 of 2-part Julian Date, TT
    ///                    tt2: Part 2 of 2-part Julian Date, TT.</returns>
    public static Tuple<int, double, double> Taitt(double tai1, double tai2) {
        double tt1 = 0, tt2 = 0;
        var status = S_Taitt(tai1, tai2, ref tt1, ref tt2);
        return Tuple.Create(status, tt1, tt2); }
    [DllImport(DllFilename, EntryPoint = "iauTaitt", CallingConvention = CallingConvention.Cdecl)]
    static extern int S_Taitt(double tai1, double tai2, ref double tt1, ref double tt2);


    /// <summary>Sofa.Taiut1(): Time scale transformation: International Atomic Time, TAI,
    /// to Universal Time, UT1.</summary>
    /// <param name="tai1">Part 1 of 2-part Julian Date, TAI</param>
    /// <param name="tai2">Part 2 of 2-part Julian Date, TAI</param>
    /// <param name="dta">UT1-TAI (seconds).</param>
    /// <returns>3-Tuple:  status [int]: 0 -> OK
    ///                    ut11: Part 1 of 2-part Julian Date, UT1
    ///                    ut12: Part 2 of 2-part Julian Date, UT1.</returns>
    public static Tuple<int, double, double> Taiut1(double tai1, double tai2, double dta) {
        double ut11 = 0, ut12 = 0;
        var status = S_Taiut1(tai1, tai2, dta, ref ut11, ref ut12);
        return Tuple.Create(status, ut11, ut12); }
    [DllImport(DllFilename, EntryPoint = "iauTaiut1", CallingConvention = CallingConvention.Cdecl)]
    static extern int S_Taiut1(double tai1, double tai2, double dta, ref double ut11, ref double ut12);


    /// <summary>Sofa.Taiutc(): Time scale transformation: International Atomic Time, TAI,
    /// to Coordinated Universal Time, UTC.</summary>
    /// <param name="tai1">Part 1 of 2-part Julian Date, TAI</param>
    /// <param name="tai2">Part 2 of 2-part Julian Date, TAI</param>
    /// <returns>3-Tuple:  status [int]: -1 -> dubious year
    ///                                   0 -> OK
    ///                                  +1 -> unacceptable date
    ///                    utc1: Part 1 of 2-part Julian Date, UTC
    ///                    utc2: Part 2 of 2-part Julian Date, UTC.</returns>
    public static Tuple<int, double, double> Taiutc(double tai1, double tai2) {
        double utc1 = 0, utc2 = 0;
        var status = S_Taiutc(tai1, tai2, ref utc1, ref utc2);
        return Tuple.Create(status, utc1, utc2); }
    [DllImport(DllFilename, EntryPoint = "iauTaiutc", CallingConvention = CallingConvention.Cdecl)]
    static extern int S_Taiutc(double tai1, double tai2, ref double utc1, ref double utc2);


    /// <summary>Sofa.Tcbtdb(): Time scale transformation: Barycentric Coordinate Time, TCB,
    /// to Barycentric Dynamical Time, TDB (~= Teph for JPL ephemerides).</summary>
    /// <param name="tcb1">Part 1 of 2-part Julian Date, TCB</param>
    /// <param name="tcb2">Part 2 of 2-part Julian Date, TCB</param>
    /// <returns>3-Tuple:  status [int]: 0 -> OK
    ///                    tdb1: Part 1 of 2-part Julian Date, TDB
    ///                    tdb2: Part 2 of 2-part Julian Date, TDB.</returns>
    public static Tuple<int, double, double> Tcbtdb(double tcb1, double tcb2) {
        double tdb1 = 0, tdb2 = 0;
        var status = S_Tcbtdb(tcb1, tcb2, ref tdb1, ref tdb2);
        return Tuple.Create(status, tdb1, tdb2); }
    [DllImport(DllFilename, EntryPoint = "iauTcbtdb", CallingConvention = CallingConvention.Cdecl)]
    static extern int S_Tcbtdb(double tcb1, double tcb2, ref double tdb1, ref double tdb2);


    /// <summary>Sofa.Tcgtt(): Time scale transformation: Geocentric Coordinate Time, TCG,
    /// to Terrestrial Time, TT.</summary>
    /// <param name="tcg1">Part 1 of 2-part Julian Date, TCG</param>
    /// <param name="tcg2">Part 2 of 2-part Julian Date, TCG</param>
    /// <returns>3-Tuple:  status [int]: 0 -> OK
    ///                    tt1: Part 1 of 2-part Julian Date, TT
    ///                    tt2: Part 2 of 2-part Julian Date, TT.</returns>
    public static Tuple<int, double, double> Tcgtt(double tcg1, double tcg2) {
        double tt1 = 0, tt2 = 0;
        var status = S_Tcgtt(tcg1, tcg2, ref tt1, ref tt2);
        return Tuple.Create(status, tt1, tt2); }
    [DllImport(DllFilename, EntryPoint = "iauTcgtt", CallingConvention = CallingConvention.Cdecl)]
    static extern int S_Tcgtt(double tcg1, double tcg2, ref double tt1, ref double tt2);


    /// <summary>Sofa.Tdbtcb(): Time scale transformation: Barycentric Dynamical Time, TDB,
    /// to Barycentric Coordinate Time, TCB.</summary>
    /// <param name="tdb1">Part 1 of 2-part Julian Date, TDB</param>
    /// <param name="tdb2">Part 2 of 2-part Julian Date, TDB</param>
    /// <returns>3-Tuple:  status [int]: 0 -> OK
    ///                    tcb1: Part 1 of 2-part Julian Date, TCB
    ///                    tcb2: Part 2 of 2-part Julian Date, TCB.</returns>
    public static Tuple<int, double, double> Tdbtcb(double tdb1, double tdb2) {
        double tcb1 = 0, tcb2 = 0;
        var status = S_Tdbtcb(tdb1, tdb2, ref tcb1, ref tcb2);
        return Tuple.Create(status, tcb1, tcb2); }
    [DllImport(DllFilename, EntryPoint = "iauTdbtcb", CallingConvention = CallingConvention.Cdecl)]
    static extern int S_Tdbtcb(double tdb1, double tdb2, ref double tcb1, ref double tcb2);


    /// <summary> Sofa.Tdbtt(): Time scale transformation: TDB (Barycentric Dynamical Time)
    /// to TT (Terrestrial Time). (TDB is essentially the same as Teph, the time argument for
    /// the JPL solar system ephemerides.)</summary>
    /// <param name="tdb1">Part 1 of 2-part Julian Date, TDB</param>
    /// <param name="tdb2">Part 2 of 2-part Julian Date, TDB</param>
    /// <param name="dtr">TDB-TT (seconds). Dominated by annual term of 1.7 milliseconds amplitude</param>
    /// <returns>3-Tuple:  status [int]: 0 -> OK
    ///                    tt1: Part 1 of 2-part Julian Date, TT
    ///                    tt2: Part 2 of 2-part Julian Date, TT.</returns>
    public static Tuple<int, double, double> Tdbtt(double tdb1, double tdb2, double dtr) {
        double tt1 = 0, tt2 = 0;
        var status = S_Tdbtt(tdb1, tdb2, dtr, ref tt1, ref tt2);
        return Tuple.Create(status, tt1, tt2); }
    [DllImport(DllFilename, EntryPoint = "iauTdbtt", CallingConvention = CallingConvention.Cdecl)]
    static extern int S_Tdbtt(double tdb1, double tdb2, double dtr, ref double tt1, ref double tt2);

    
    /// <summary>Sofa.Tttai(): Time scale transformation: Terrestrial Time, TT,
    /// to International Atomic Time, TAI.</summary>
    /// <param name="tt1">Part 1 of 2-part Julian Date, TT</param>
    /// <param name="tt2">Part 2 of 2-part Julian Date, TT</param>
    /// <returns>3-Tuple:  status [int]: 0 -> OK
    ///                    tai1: Part 1 of 2-part Julian Date, TAI
    ///                    tai2: Part 2 of 2-part Julian Date, TAI.</returns>
    public static Tuple<int, double, double> Tttai(double tt1, double tt2) {
        double tai1 = 0, tai2 = 0;
        var status = S_Tttai(tt1, tt2, ref tai1, ref tai2);
        return Tuple.Create(status, tai1, tai2); }
    [DllImport(DllFilename, EntryPoint = "iauTttai", CallingConvention = CallingConvention.Cdecl)]
    static extern int S_Tttai(double tt1, double tt2, ref double tai1, ref double tai2);


    /// <summary>Sofa.Tttcg(): Time scale transformation: Terrestrial Time, TT,
    /// to Geocentric Coordinate Time, TCG.</summary>
    /// <param name="tt1">Part 1 of 2-part Julian Date, TT</param>
    /// <param name="tt2">Part 2 of 2-part Julian Date, TT</param>
    /// <returns>3-Tuple:  status [int]: 0 -> OK
    ///                    tcg1: Part 1 of 2-part Julian Date, TCG
    ///                    tcg2: Part 2 of 2-part Julian Date, TCG</returns>
    public static Tuple<int, double, double> Tttcg(double tt1, double tt2) {
        double tcg1 = 0, tcg2 = 0;
        var status = S_Tttcg(tt1, tt2, ref tcg1, ref tcg2);
        return Tuple.Create(status, tcg1, tcg2); }
    [DllImport(DllFilename, EntryPoint = "iauTttcg", CallingConvention = CallingConvention.Cdecl)]
    static extern int S_Tttcg(double tt1, double tt2, ref double tcg1, ref double tcg2);


    /// <summary> Sofa.Tttdb(): Time scale transformation: Terrestrial Time, TT,
    /// to Barycentric Dynamical Time, TDB.</summary>
    /// <param name="tt1">Part 1 of 2-part Julian Date, TT</param>
    /// <param name="tt2">Part 2 of 2-part Julian Date, TT</param>
    /// <param name="dtr">TDB-TT (seconds). Dominated by annual term of 1.7 milliseconds amplitude.</param>
    /// <returns>3-Tuple:  status [int]: 0 -> OK
    ///                    tdb1: Part 1 of 2-part Julian Date, TDB
    ///                    tdb2: Part 2 of 2-part Julian Date, TDB</returns>
    public static Tuple<int, double, double> Tttdb(double tt1, double tt2, double dtr) {
        double tdb1 = 0, tdb2 = 0;
        var status = S_Tttdb(tt1, tt2, dtr, ref tdb1, ref tdb2);
        return Tuple.Create(status, tdb1, tdb2); }
    [DllImport(DllFilename, EntryPoint = "iauTttdb", CallingConvention = CallingConvention.Cdecl)]
    static extern int S_Tttdb(double tt1, double tt2, double dtr, ref double tdb1, ref double tdb2);


    /// <summary>Sofa.Ttut1(): Time scale transformation: Terrestrial Time, TT,
    /// to Universal Time, UT1.</summary>
    /// <param name="tt1">Part 1 of 2-part Julian Date, TT</param>
    /// <param name="tt2">Part 2 of 2-part Julian Date, TT</param>
    /// <param name="dt">TT-UT1 (seconds).</param>
    /// <returns>3-Tuple:  status [int]: 0 -> OK
    ///                    ut11: Part 1 of 2-part Julian Date, UT1
    ///                    ut12: Part 2 of 2-part Julian Date, UT1.</returns>
    public static Tuple<int, double, double> Ttut1(double tt1, double tt2, double dt) {
        double ut11 = 0, ut12 = 0;
        var status = S_Ttut1(tt1, tt2, dt, ref ut11, ref ut12);
        return Tuple.Create(status, ut11, ut12); }
    [DllImport(DllFilename, EntryPoint = "iauTtut1", CallingConvention = CallingConvention.Cdecl)]
    static extern int S_Ttut1(double tt1, double tt2, double dt, ref double ut11, ref double ut12);


    /// <summary>Sofa.Ut1tai(): TTime scale transformation: Universal Time, UT1,
    /// to International Atomic Time, TAI.</summary>
    /// <param name="ut11">Part 1 of 2-part Julian Date, UT1</param>
    /// <param name="ut12">Part 2 of 2-part Julian Date, UT1</param>
    /// <param name="dta">UT1-TAI, available IERS tabulations (seconds).</param>
    /// <returns>3-Tuple:  status [int]: 0 -> OK
    ///                    tai1: Part 1 of 2-part Julian Date, TAI
    ///                    tai2: Part 2 of 2-part Julian Date, TAI.</returns>
    public static Tuple<int, double, double> Ut1tai(double ut11, double ut12, double dta) {
        double tai1 = 0, tai2 = 0;
        var status = S_Ut1tai(ut11, ut12, dta, ref tai1, ref tai2);
        return Tuple.Create(status, tai1, tai2); }
    [DllImport(DllFilename, EntryPoint = "iauUt1tai", CallingConvention = CallingConvention.Cdecl)]
    static extern int S_Ut1tai(double ut11, double ut12, double dta, ref double tai1, ref double tai2);


    /// <summary>Sofa.Ut1tt(): Time scale transformation: Universal Time, UT1,
    /// to Terrestrial Time, TT.</summary>
    /// <param name="ut11">Part 1 of 2-part Julian Date, UT1</param>
    /// <param name="ut12">Part 2 of 2-part Julian Date, UT1</param>
    /// <param name="dt">TT-UT1, "classical Delta T" (seconds)</param>
    /// <returns>3-Tuple:  status [int]: 0 -> OK
    ///                    tt1: Part 1 of 2-part Julian Date, TT
    ///                    tt2: Part 2 of 2-part Julian Date, TT.</returns>
    public static Tuple<int, double, double> Ut1tt(double ut11, double ut12, double dt) {
        double tt1 = 0, tt2 = 0;
        var status = S_Ut1tt(ut11, ut12, dt, ref tt1, ref tt2);
        return Tuple.Create(status, tt1, tt2); }
    [DllImport(DllFilename, EntryPoint = "iauUt1tt", CallingConvention = CallingConvention.Cdecl)]
    static extern int S_Ut1tt(double ut11, double ut12, double dt, ref double tt1, ref double tt2);


    /// <summary>Sofa.Ut1utc(): Time scale transformation: Universal Time, UT1,
    /// to Coordinated Universal Time, UTC.</summary>
    /// <param name="ut11">Part 1 of 2-part Julian Date, UT1</param>
    /// <param name="ut12">Part 2 of 2-part Julian Date, UT1</param>
    /// <param name="dut1">UT1-UTC, "Delta UT1", available IERS tabulations (seconds)</param>
    /// <returns>3-Tuple:  status [int]: -1 -> dubious year, 0 -> OK, +1 -> unacceptable date
    ///                    utc1: Part 1 of 2-part Julian Date, UTC
    ///                    utc2: Part 2 of 2-part Julian Date, UTC.</returns>
    public static Tuple<int, double, double> Ut1utc(double ut11, double ut12, double dut1) {
        double utc1 = 0, utc2 = 0;
        var status = S_Ut1utc(ut11, ut12, dut1, ref utc1, ref utc2);
        return Tuple.Create(status, utc1, utc2); }
    [DllImport(DllFilename, EntryPoint = "iauUt1utc", CallingConvention = CallingConvention.Cdecl)]
    static extern int S_Ut1utc(double ut11, double ut12, double dut1, ref double utc1, ref double utc2);


    /// <summary>Sofa.Utctai(): Time scale transformation: Coordinated Universal Time, UTC,
    /// to International Atomic Time, TAI.</summary>
    /// <param name="utc1">Part 1 of 2-part Julian Date, UTC</param>
    /// <param name="utc2">Part 2 of 2-part Julian Date, UTC</param>
    /// <returns>3-Tuple:  status [int]: +1 -> dubious year, 0 -> OK, -1 -> unacceptable date
    ///                    tai1: Part 1 of 2-part Julian Date, TAI
    ///                    tai2: Part 2 of 2-part Julian Date, TAI.</returns>
    public static Tuple<int, double, double> Utctai(double utc1, double utc2) {
        double tai1 = 0, tai2 = 0;
        var status = S_Utctai(utc1, utc2, ref tai1, ref tai2);
        return Tuple.Create(status, tai1, tai2); }
    [DllImport(DllFilename, EntryPoint = "iauUtctai", CallingConvention = CallingConvention.Cdecl)]
    static extern int S_Utctai(double utc1, double utc2, ref double tai1, ref double tai2);

    
    /// <summary>Sofa.Utcut1(): Time scale transformation: Coordinated Universal Time, UTC,
    /// to Universal Time, UT1.</summary>
    /// <param name="utc1">Part 1 of 2-part Julian Date, UTC</param>
    /// <param name="utc2">Part 2 of 2-part Julian Date, UTC</param>
    /// <param name="dut1">UT1-UTC, "Delta UT1", available IERS tabulations (seconds)</param>
    /// <returns>3-Tuple:  status [int]: +1 -> dubious year, 0 -> OK, -1 -> unacceptable date
    ///                    ut11: Part 1 of 2-part Julian Date, UT1
    ///                    ut12: Part 2 of 2-part Julian Date, UT1.</returns>
    public static Tuple<int, double, double> Utcut1(double utc1, double utc2, double dut1) {
        double ut11 = 0, ut12 = 0;
        var status = S_Utcut1(utc1, utc2, dut1, ref ut11, ref ut12);
        return Tuple.Create(status, ut11, ut12); }
    [DllImport(DllFilename, EntryPoint = "iauUtcut1", CallingConvention = CallingConvention.Cdecl)]
    static extern int S_Utcut1(double utc1, double utc2, double dut1, ref double ut11, ref double ut12);

    #endregion  //  Time Scales.


    #region "Astronomy: Ecliptic Coordinates"
    /*######### Astronomy: Ecliptic Coordinates ##########################################################*/
    
    /// <summary>Sofa.Eceq06(): Transformation from ecliptic coordinates (mean equinox and ecliptic
    /// of date) to ICRS RA,Dec, using the IAU 2006 precession model.</summary>
    /// <param name="date1">Part 1 of 2-Part TT quasi-Julian date</param>
    /// <param name="date2">Part 2 of 2-Part TT quasi-Julian date</param>
    /// <param name="dl">Ecliptic longitude (radians)</param>
    /// <param name="db">Ecliptic latitude (radians)</param>
    /// <returns>2-tuple of doubles: dr: ICRS Right Ascension (radians)
    ///                              dd: ICRS Declination (radians)</returns>
    public static Tuple<double, double> Eceq06(double date1, double date2, double dl, double db) {
        double dr = 0, dd = 0;
        S_Eceq06(date1, date2, dl, db, ref dr, ref dd);
        return Tuple.Create(dr, dd); }
    [DllImport(DllFilename, EntryPoint = "iauEceq06", CallingConvention = CallingConvention.Cdecl)]
    static extern void S_Eceq06(double date1, double date2, double dl, double db,
        ref double dr, ref double dd);


    /// <summary> Sofa.Ecm06(): Generate rotation matrix, ICRS equatorial to Ecliptic, IAU 2006.</summary>
    /// <param name="date1">TT as 2-part Julian date, part 1</param>
    /// <param name="date2">TT as 2-part Julian date, part 2</param>
    /// <returns>Rotation matrix, ICRS to Ecliptic (3x3 matrix of doubles)</returns>
    public static double[,] Ecm06(double date1, double date2) {
        var rm = new double[3, 3] {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
        S_Ecm06(date1, date2, rm);
        return rm; }
    [DllImport(DllFilename, EntryPoint = "iauEcm06", CallingConvention = CallingConvention.Cdecl)]
    static extern void S_Ecm06(double date1, double date2, [In, Out] double[,] rm);


    /// <summary>Sofa.Eqec06(): Transformation from ICRS equatorial coordinates to ecliptic coordinates
    /// (mean equinox and ecliptic of date), using the IAU 2006 precession model. </summary>
    /// <param name="date1">Part 1 of 2-Part TT quasi-Julian date</param>
    /// <param name="date2">Part 2 of 2-Part TT quasi-Julian date</param>
    /// <param name="dr">ICRS Right Ascension (radians)</param>
    /// <param name="dd">ICRS Declination (radians)</param>
    /// <returns>2-tuple of doubles: dl: Ecliptic Longitude (radians)
    ///                              db: Ecliptic latitude (radians)</returns>
    public static Tuple<double, double> Eqec06(double date1, double date2, double dr, double dd) {
        double dl = 0, db = 0;
        S_Eqec06(date1, date2, dr, dd, ref dl, ref db);
        return Tuple.Create(dl, db); }
    [DllImport(DllFilename, EntryPoint = "iauEqec06", CallingConvention = CallingConvention.Cdecl)]
    static extern void S_Eqec06(double date1, double date2, double dr, double dd,
        ref double dl, ref double db);
    
    // Lteceq():  NOT IMPLEMENTED: Very long-term coordinates not relevant to amateur astronomy.
    // Ltecm():   NOT IMPLEMENTED: Very long-term coordinates not relevant to amateur astronomy.
    // Lteqec():  NOT IMPLEMENTED: Very long-term coordinates not relevant to amateur astronomy.

    #endregion //  Ecliptic Coordinates.


    #region "Astronomy: Earth Rotation and Sidereal Time"
    /*######### Astronomy: Earth Rotation and Sidereal Time ##############################################*/

    // Ee00():      NOT IMPLEMENTED: AstroLib favors CIO over equation of equinoxes.
    // Ee00a():     NOT IMPLEMENTED: AstroLib favors CIO over equation of equinoxes.
    // Ee00b():     NOT IMPLEMENTED: AstroLib favors CIO over equation of equinoxes.
    // Ee06a():     NOT IMPLEMENTED: AstroLib favors CIO over equation of equinoxes.
    // Eect00():    NOT IMPLEMENTED: AstroLib favors CIO over equation of equinoxes.
    // Eqeq94():    NOT IMPLEMENTED: AstroLib favors CIO over equation of equinoxes.
    
    /// <summary>Sofa.Era00(): Earth rotation angle (IAU 2000 model).</summary>
    /// <param name="date1">Part 1 of Julian Date, UTC</param>
    /// <param name="date2">Part 2 of Julian Date, UTC</param>
    /// <returns>Earth Rotation Angle (radians, range [0, 2*pi])</returns>.
    public static double Era00(double date1, double date2) {
        var era = S_Era00(date1, date2);
        return era; }
    [DllImport(DllFilename, EntryPoint = "iauEra00", CallingConvention = CallingConvention.Cdecl)]
    static extern double S_Era00(double date1, double date2);

    // Gmst00():    NOT IMPLEMENTED: Prefer Gmst06() with its newer model.
    
    /// <summary>Sofa.Gmst06(): Greenwich mean sidereal time (consistent with IAU 2006 precession).
    /// Per SOFA manual "...if UT1 is used for both purposes, errors of order 100 microarcseconds result",
    ///     and such errors are almost always acceptable, we preserve herein the option to provide both.
    /// </summary>
    /// <param name="uta">Part 1 of Julian Date, UT1</param>
    /// <param name="utb">Part 2 of Julian Date, UT1</param>
    /// <param name="tta">Part 1 of Julian Date, TT (or use UT1 for milliarcsecond accuracy)</param>
    /// <param name="ttb">Part 2 of Julian Date, TT (or use UT1 for milliarcsecond accuracy)</param>
    /// <returns>Greenwich mean sidereal time (radians)</returns>.
    public static double Gmst06(double uta, double utb, double tta, double ttb) {
        var gmst = S_Gmst06(uta, utb, tta, ttb);
        return gmst; }
    [DllImport(DllFilename, EntryPoint = "iauGmst06", CallingConvention = CallingConvention.Cdecl)]
    static extern double S_Gmst06(double uta, double utb, double tta, double ttb);
    
    // Gmst82():    NOT IMPLEMENTED: Prefer Gmst06() with its newer model.
    // Gst00a():    NOT IMPLEMENTED: Prefer Gst06a() with its newer model.
    // Gst00b():    NOT IMPLEMENTED: Prefer Gst06a() with its newer model.
    // Gst06():     NOT IMPLEMENTED: Prefer Gst06a() which does not require NPB matrix.
    
    /// <summary>Sofa.Gst06a(): Greenwich apparent sidereal time
    /// (consistent with IAU 2000 and 2006 resolutions).</summary>
    /// <param name="uta">Part 1 of Julian Date, UT1</param>
    /// <param name="utb">Part 2 of Julian Date, UT1</param>
    /// <param name="tta">Part 1 of Julian Date, TT (or use UT1 for milliarcsecond accuracy)</param>
    /// <param name="ttb">Part 2 of Julian Date, TT (or use UT1 for milliarcsecond accuracy)</param>
    /// <returns>Greenwich mean sidereal time (radians)</returns>.
    public static double Gst06a(double uta, double utb, double tta, double ttb) {
        var gast = S_Gst06a(uta, utb, tta, ttb);
        return gast; }
    [DllImport(DllFilename, EntryPoint = "iauGst06a", CallingConvention = CallingConvention.Cdecl)]
    static extern double S_Gst06a(double uta, double utb, double tta, double ttb);

    // Gst94():    NOT IMPLEMENTED: Prefer Gst06a() with its newer model.

    #endregion  //  Earth Rotation and Sidereal Time.


    #region "Astronomy: Ephemerides"
    /*######### Astronomy: Ephemerides ###################################################################*/
    
    /// <summary>Sofa.Epv00(): Earth position and velocity, both heliocentric and barycentric,
    /// with respect to the Barycentric Celestial Reference System.</summary>
    /// <param name="date1">Part 1 of quasi-Julian form of TDB (or possibly TT) date</param>
    /// <param name="date2">Part 2 of quasi-Julian form of TDB (or possibly TT) date</param>
    /// <returns>3-tuple:
    ///         status [int]: 0 -> OK, +1 -> date outside [1900-2100 AD] 
    ///         pvh [double[2][3]]: heliocentric Earth position and velocity (AU, AU/day)
    ///         pvb [double[2][3]]: barycentric Earth position and velocity (AU, AU/day) </returns>.
    public static Tuple<int, double[,], double[,]> Epv00(double date1, double date2) {
        var pvh = new double[2, 3] {{0, 0, 0}, {0, 0, 0}};
        var pvb = new double[2, 3] {{0, 0, 0}, {0, 0, 0}};
        var status = S_Epv00(date1, date2, pvh, pvb);
        return Tuple.Create(status, pvh, pvb); }
    [DllImport(DllFilename, EntryPoint = "iauEpv00", CallingConvention = CallingConvention.Cdecl)]
    static extern int S_Epv00(double date1, double date2, [In, Out] double[,] pvh, [In, Out] double[,] pvb);


    /// <summary>Sofa.Moon98(): Approximate geocentric position and velocity of the Moon.</summary>
    /// <param name="date1">Part 1 of quasi-Julian form of TT (or possibly TDB) date</param>
    /// <param name="date2">Part 2 of quasi-Julian form of TT (or possibly TDB) date</param>
    /// <returns> pv [double[2][3]]: Moon position and velocity, GCRS (in AU, AU/d)</returns>
    public static double[,] Moon98(double date1, double date2) {
        var pv = new double[2, 3] {{0, 0, 0}, {0, 0, 0}};
        S_Moon98(date1, date2, pv);
        return pv; }
    [DllImport(DllFilename, EntryPoint = "iauMoon98", CallingConvention = CallingConvention.Cdecl)]
    static extern void S_Moon98(double date1, double date2, [In, Out] double[,] pv);


    /// <summary>Sofa.Plan94(): Approximate heliocentric position and velocity of a nominated major
    /// planet: Mercury, Venus, EMB, Mars, Jupiter, Saturn, Uranus or
    /// Neptune (but not the Earth itself).</summary>
    /// <param name="date1">Part 1 of quasi-Julian form of TDB (or possibly TT) date</param>
    /// <param name="date2">Part 2 of quasi-Julian form of TDB (or possibly TT) date</param>
    /// <param name="np">planet (1=Mercury, 2=Venus, 3=EMB (Earth-Moon Barycenter; for earth alone use
    /// Epv00), 4=Mars, 5=Jupiter, 6=Saturn, 7=Uranus, 8=Neptune)</param>
    /// <returns>2-tuple: status [int]: -1 -> error: illegal planet number np (must be 1-8)
    ///                                  0 -> OK
    ///                                 +1 -> warning: year outside 1000-3000
    ///                                 +2 -> warning: failed to converge
    ///                   pv [double[2][3]]: Planet position and velocity,
    ///                                 heliocentric J2000.0 (in AU, AU/d)</returns>.
    public static Tuple<int, double[,]> Plan94(double date1, double date2, int np) {
        var pv = new double[2, 3] {{0, 0, 0}, {0, 0, 0}};
        var status = S_Plan94(date1, date2, np, pv);
        return Tuple.Create(status, pv); }
    [DllImport(DllFilename, EntryPoint = "iauPlan94", CallingConvention = CallingConvention.Cdecl)]
    static extern int S_Plan94(double date1, double date2, int np, [In, Out] double[,] pv);

    #endregion  //  Ephemerides.


    #region "Astronomy: Fundamental Arguments"
    /*######### Astronomy: Fundamental Arguments #########################################################*/

    // All 14 functions in Fundamental Arguments section:  NOT IMPLEMENTED.

    #endregion  //  Fundamental Arguments.


    #region "Astronomy: Galactic Coordinates"
    /*######### Astronomy: Galactic Coordinates ##########################################################*/
    
    /// <summary>Sofa.G2icrs(): Transformation from Galactic Coordinates (1958) to ICRS.</summary>
    /// <param name="dl">Galactic longitude (radians)</param>
    /// <param name="db">Galactic latitude (radians)</param>
    /// <returns>2-tuple of doubles: dr: ICRS right ascension (radians)
    ///                              dd: ICRS declination (radians) </returns>.
    public static Tuple<double, double> G2icrs(double dl, double db) {
        double dr = 0, dd = 0;
        S_G2icrs(dl, db, ref dr, ref dd);
        return Tuple.Create(dr, dd); }
    [DllImport(DllFilename, EntryPoint = "iauG2icrs", CallingConvention = CallingConvention.Cdecl)]
    static extern void S_G2icrs(double dl, double db, ref double dr, ref double dd);


    /// <summary>Sofa.Icrs2g(): Transformation from ICRS to Galactic Coordinates (1958).</summary>
    /// <param name="dr">ICRS right ascension (radians)</param>
    /// <param name="dd">ICRS Declination (radians)</param>
    /// <returns>2-tuple of doubles: dl: Galactic longitude (radians)
    ///                              db: Galactic latitude (radians) </returns>.
    public static Tuple<double, double> Icrs2g(double dr, double dd) {
        double dl = 0, db = 0;
        S_Icrs2g(dr, dd, ref dl, ref db);
        return Tuple.Create(dl, db); }
    [DllImport(DllFilename, EntryPoint = "iauIcrs2g", CallingConvention = CallingConvention.Cdecl)]
    static extern void S_Icrs2g(double dr, double dd, ref double dl, ref double db);

    #endregion  //  Galactic Coordinates.


    #region "Astronomy: Geocentric/Geodetic Transformations"
    /*######### Astronomy: Geocentric/Geodetic Transformations ###########################################*/

    // Eform():  NOT IMPLEMENTED:  No need for customized Earth reference ellipsoids; AstroLib uses WGS84.

    /// <summary>Sofa.Gc2gd(): Transform geocentric coordinates to geodetic using the specified
    /// reference ellipsoid. NB: almost always use n=1 for WGS84.</summary>
    /// <param name="n">Ellipsoid identifier; n=1 -> WGS84, the normal case</param>
    /// <param name="xyz">Geocentric vector (all 3 values in meters)</param>
    /// <returns>4-tuple:  status [int]: 0 -> OK, -1 -> illegal identifier (n), -2 -> internal error
    ///                    elong: Longitude (radians)
    ///                    phi: Latitude (radians)
    ///                    height: Height above ellipsoid (meters, geodetic). </returns>
    public static Tuple<int, double, double, double> Gc2gd(int n, double[] xyz) {
        double elong = 0, phi = 0, height = 0;
        var status = S_Gc2gd(n, xyz, ref elong, ref phi, ref height);
        return Tuple.Create(status, elong, phi, height); }
    [DllImport(DllFilename, EntryPoint = "iauGc2gd", CallingConvention = CallingConvention.Cdecl)]
    static extern int S_Gc2gd(int n, double[] xyz, ref double elong, ref double phi, ref double height);

    // Gc2gde():  NOT IMPLEMENTED:  No need for customized Earth reference ellipsoids; AstroLib uses WGS84.

    /// <summary>Sofa.Gd2gc(): Transform geodetic coordinates to geocentric using the specified
    /// reference ellipsoid.</summary>
    /// <param name="n">Ellipsoid identifier; n=1 -> WGS84, the normal case.</param>
    /// <param name="elong">Longitude (radians).</param>
    /// <param name="phi">Latitude (radians).</param>
    /// <param name="height">Height above ellipsoid (meters, geodetic).</param>
    /// <returns>2-Tuple:  status [int]: 0 -> OK, -1 -> illegal identifier (n), -2 -> illegal case
    ///                    xyz [double[3]]: Geocentric vector (meters).</returns>
    public static Tuple<int, double[]> Gd2gc(int n, double elong, double phi, double height) {
        var xyz = new double[3] {0, 0, 0};
        var status = S_Gd2gc(n, elong, phi, height, xyz);
        return Tuple.Create(status, xyz); }
    [DllImport(DllFilename, EntryPoint = "iauGd2gc", CallingConvention = CallingConvention.Cdecl)]
    static extern int S_Gd2gc(int n, double elong, double phi, double height, [In, Out] double[] xyz);

    // Gd2gce():  NOT IMPLEMENTED:  No need for customized Earth reference ellipsoids; AstroLib uses WGS84.

    #endregion  //  Geocentric/Geodetic Transformations.


    #region "Astronomy: Gnomonic Projections"
    /*######### Astronomy: Gnomonic Projections ##########################################################*/
    
    // Tpors():  NOT IMPLEMENTED: Finding tangent point not a common need in amateur astronomy.
    // Tporv():  NOT IMPLEMENTED: Finding tangent point not a common need in amateur astronomy.

    /// <summary> Sofa.Tpsts(): Project a tangent plane to a star's celestial spherical coordinates.
    /// That is, in the tangent plane projection, given the star’s rectangular coordinates
    /// and the spherical coordinates of the tangent point, solve for the
    /// spherical coordinates of the star.</summary>
    /// <param name="xi">Rectangular coordinates of star image (radians at tangent point)</param>
    /// <param name="eta">Rectangular coordinates of star image, due North (radians @ tangent pt)</param>
    /// <param name="a0">Tangent point's spherical coordinates (radians)</param>
    /// <param name="b0">Tangent point's spherical coordinates (radians)</param>
    /// <returns>2-tuple of doubles: a: Star's spherical coordinates (radians).
    ///                              b: Star's spherical coordinates (radians).</returns>.
    public static Tuple<double, double> Tpsts(double xi, double eta, double a0, double b0) {
        double a = 0, b = 0;
        S_Tpsts(xi, eta, a0, b0, ref a, ref b);
        return Tuple.Create(a, b); }
    [DllImport(DllFilename, EntryPoint = "iauTpsts", CallingConvention = CallingConvention.Cdecl)]
    static extern void S_Tpsts(double xi, double eta, double a0, double b0, ref double a, ref double b);

    
    /// <summary> Sofa.Tpstv(): Project a tangent plane to a star's celestial unit vector.
    /// That is, in the tangent plane projection, given the star’s rectangular coordinates
    /// and the unit vector (direction cosines) of the tangent point, solve for the
    /// unit vector (direction cosines) of the star.</summary>
    /// <param name="xi">Rectangular coordinates of star image (radians at tangent point)</param>
    /// <param name="eta">Rectangular coordinates of star image, due North (radians @ tangent pt)</param>
    /// <param name="v0">Tangent point's unit vector (direction cosines)</param>
    /// <returns>Star's unit vector (direction cosines)</returns>
    public static double[] Tpstv(double xi, double eta, double[] v0) {
        var v  = new double[3] {0, 0, 0};
        S_Tpstv(xi, eta, v0, v);
        return v; }
    [DllImport(DllFilename, EntryPoint = "iauTpstv", CallingConvention = CallingConvention.Cdecl)]
    static extern void S_Tpstv (double xi, double eta, double[] v0, [In, Out] double[] v);


    /// <summary>Sofa.Tpxes(): Project a star's celestial spherical coordinates to find its rectangular
    /// coordinates. That is, in the tangent plane projection, given celestial spherical coordinates for
    /// a star and the tangent point, solve for the star’s rectangular coordinates
    /// in the tangent plane.</summary>
    /// <param name="a">Star's spherical coordinates (radians)</param>
    /// <param name="b">Star's spherical coordinates (radians)</param>
    /// <param name="a0">Tangent point's spherical coordinates (radians)</param>
    /// <param name="b0">Tangent point's spherical coordinates (radians)</param>
    /// <returns>3-tuple:
    ///          status [int]: 0 -> OK
    ///                        1 -> star too far from axis
    ///                        2 -> antistar on tangent plane
    ///                        3 -> antistar too far from axis
    ///          xi:  Rectangular coordinates of star image (radians at tangent point).
    ///          eta: Rectangular coordinates of star image, due North (radians at tangent pt).</returns>.
    public static Tuple<int, double, double> Tpxes(double a, double b, double a0, double b0) {
        double xi = 0, eta = 0;
        var status = S_Tpxes(a, b, a0, b0, ref xi, ref eta);
        return Tuple.Create(status, xi, eta); }
    [DllImport(DllFilename, EntryPoint = "iauTpxes", CallingConvention = CallingConvention.Cdecl)]
    static extern int S_Tpxes(double a, double b, double a0, double b0, ref double xi, ref double eta);

    
    /// <summary> Sofa.Tpxev():  Project a star's unit vector (direction cosines) to find its
    /// rectangular coordinates. That is, in the tangent plane projection,
    /// given a star's celestial unit vector (direction cosines) and the tangent point's
    /// unit vector (direction cosines), solve for the star's rectangular coordinates
    /// in the tangent plane.</summary>
    /// <param name="v">Star's unit vector (direction cosines).</param>
    /// <param name="v0">Tangent point's unit vector (direction cosines).</param>
    /// <returns>3-tuple:
    ///          status [int]: 0 -> OK
    ///                        1 -> star too far from axis
    ///                        2 -> antistar on tangent plane
    ///                        3 -> antistar too far from axis
    ///          xi:  Rectangular coordinates of star image (radians at tangent point).
    ///          eta: Rectangular coordinates of star image, due North (radians at tangent pt).</returns>.
    public static Tuple<int, double, double> Tpxev(double[] v, double[] v0) {
        double xi = 0, eta = 0;
        var status = S_Tpxev(v, v0, ref xi, ref eta);
        return Tuple.Create(status, xi, eta); }
    [DllImport(DllFilename, EntryPoint = "iauTpxev", CallingConvention = CallingConvention.Cdecl)]
    static extern int S_Tpxev(double[] v, double[] v0, ref double xi, ref double eta);
    
    #endregion  //  Gnomonic Projections.


    #region "Astronomy: Horizon/Equatorial Coordinates"
    /*######### Astronomy: Horizon/Equatorial Coordinates ################################################*/
    
    /// <summary>Sofa.Ae2hd(): Transform (azimuth, altitude) to (hour angle, declination).</summary>
    /// <param name="az">Azimuth (radians).</param>
    /// <param name="el">Altitude (elevation, radians).</param>
    /// <param name="phi">Site latitude (radians).</param>
    /// <returns>2-tuple of doubles: ha: Hour angle at site (radians)
    ///                              dec: Declination at site (radians)</returns>.
    public static Tuple<double, double> Ae2hd(double az, double el, double phi) {
        double ha = 0, dec = 0;
        S_Ae2hd(az, el, phi, ref ha, ref dec);
        return Tuple.Create(ha, dec); }
    [DllImport(DllFilename, EntryPoint = "iauAe2hd", CallingConvention = CallingConvention.Cdecl)]
    static extern void S_Ae2hd(double az, double el, double phi, ref double ha, ref double dec);


    /// <summary>Sofa.Hd2ae(): Transform (hour angle, declination) to (azimuth, altitude).</summary>
    /// <param name="ha">Local hour angle (radians).</param>
    /// <param name="dec">Site declination (radians).</param>
    /// <param name="phi">Site latitude (radians).</param>
    /// <returns>2-tuple of doubles: az: Azimuth (radians)
    ///                              el: Altitude (elevation, radians)</returns>.
    public static Tuple<double, double> Hd2ae(double ha, double dec, double phi) {
        double az = 0, el = 0;
        S_Hd2ae(ha, dec, phi, ref az, ref el);
        return Tuple.Create(az, el); }
    [DllImport(DllFilename, EntryPoint = "iauHd2ae", CallingConvention = CallingConvention.Cdecl)]
    static extern void S_Hd2ae(double ha, double dec, double phi, ref double az, ref double el);
    
    // Hd2pa():  NOT IMPLEMENTED: Parallactic angle of low priority for amateur astronomy. (add later?)

    #endregion  //  Horizon/Equatorial Coordinates.


    #region "Astronomy: Precession/Nutation/Polar Motion"
    /*######### Astronomy: Precession/Nutation/Polar Motion ##############################################*/

    // Bi00():    NOT IMPLEMENTED: Frame bias not relevant to amateur astronomy. Outdated model.
    // Bp00():    NOT IMPLEMENTED: Prefer NPB matrix usage to individual components.
    // Bp06():    NOT IMPLEMENTED: Prefer NPB matrix usage to individual components.
    // Bpn2xy():  NOT IMPLEMENTED: CIP XY probably not relevant to amateur astronomy. (implement later?)
    // C2i00a():  NOT IMPLEMENTED: Prefer transforms to direct matrix usage. Outdated model.
    // C2i00b():  NOT IMPLEMENTED: Prefer transforms to direct matrix usage. Outdated model.
    // C2i06a():  NOT IMPLEMENTED: Prefer transforms to direct matrix usage.
    // C2ibpn():  NOT IMPLEMENTED: Prefer transforms to direct matrix usage.
    // C2ixy():   NOT IMPLEMENTED: Prefer transforms to direct matrix usage.
    // C2ixys():  NOT IMPLEMENTED: Prefer transforms to direct matrix usage.
    // C2t00a():  NOT IMPLEMENTED: Prefer transforms to direct matrix usage. Outdated model.
    // C2t00b():  NOT IMPLEMENTED: Prefer transforms to direct matrix usage. Outdated model.
    // C2t06a():  NOT IMPLEMENTED: Prefer transforms to direct matrix usage.
    // C2tcio():  NOT IMPLEMENTED: Prefer transforms to direct matrix usage.
    // C2teqx():  NOT IMPLEMENTED: Prefer transforms to direct matrix usage.
    // C2tpe():   NOT IMPLEMENTED: Prefer transforms to direct matrix usage.
    // C2txy():   NOT IMPLEMENTED: Prefer transforms to direct matrix usage.
    // Eo06a():   NOT IMPLEMENTED: Equation of the origins of low relevance to amateur astronomy. (later?)
    // Eors():    NOT IMPLEMENTED: Equation of the origins of low relevance to amateur astronomy. (later?)
    // Fw2m():    NOT IMPLEMENTED: Fukushima-Williams angles of low relevance to amateur astronomy.
    // Fw2xy():   NOT IMPLEMENTED: Fukushima-Williams angles of low relevance to amateur astronomy.
    // Ltp():     NOT IMPLEMENTED: Long-term precession (40K year) of low relevance to amateur astronomy.
    // Ltpb():    NOT IMPLEMENTED: Long-term precession (40K year) of low relevance to amateur astronomy.
    // Ltpecl():  NOT IMPLEMENTED: Long-term precession (40K year) of low relevance to amateur astronomy.
    // Ltpequ():  NOT IMPLEMENTED: Long-term precession (40K year) of low relevance to amateur astronomy.
    // Num00a():  NOT IMPLEMENTED: Nutation matrix of low relevance to amateur astronomy.
    // Num00b():  NOT IMPLEMENTED: Nutation matrix of low relevance to amateur astronomy.
    // Num06a():  NOT IMPLEMENTED: Nutation matrix of low relevance to amateur astronomy.
    // Numat():   NOT IMPLEMENTED: Nutation matrix of low relevance to amateur astronomy.
    // Nut00a():  NOT IMPLEMENTED: Nutation of low relevance to amateur astronomy.
    // Nut00b():  NOT IMPLEMENTED: Nutation of low relevance to amateur astronomy.
    // Nut06a():  NOT IMPLEMENTED: Nutation of low relevance to amateur astronomy.
    // Nut80():   NOT IMPLEMENTED: Nutation of low relevance to amateur astronomy. Outdated model.
    // Nutm80():  NOT IMPLEMENTED: Nutation of low relevance to amateur astronomy. Outdated model.
    // Obl06():   NOT IMPLEMENTED: Mean obliquity of low relevance to amateur astronomy.
    // Obl80():   NOT IMPLEMENTED: Mean obliquity of low relevance to amateur astronomy. Outdated model.
    // Pb06():    NOT IMPLEMENTED: Precession angles of uncertain relevance to amateur astronomy.
    // Pfw06():   NOT IMPLEMENTED: Precession angles of uncertain relevance to amateur astronomy.
    // Pmat00():  NOT IMPLEMENTED: Prefer transforms to direct matrix usage.
    // Pmat06():  NOT IMPLEMENTED: Prefer transforms to direct matrix usage.
    // Pmat76():  NOT IMPLEMENTED: Prefer transforms to direct matrix usage. Outdated model.
    // Pn00():    NOT IMPLEMENTED: Prefer transforms to direct precession-nutation usage.
    // Pn00a():   NOT IMPLEMENTED: Prefer transforms to direct precession-nutation usage.
    // Pn00b():   NOT IMPLEMENTED: Prefer transforms to direct precession-nutation usage.
    // Pn06():    NOT IMPLEMENTED: Prefer transforms to direct precession-nutation usage.
    // Pn06a():   NOT IMPLEMENTED: Prefer transforms to direct precession-nutation usage.
    // Pnm00a():  NOT IMPLEMENTED: Prefer transforms to direct matrix usage.
    // Pnm00b():  NOT IMPLEMENTED: Prefer transforms to direct matrix usage.
    // Pnm06a():  NOT IMPLEMENTED: Prefer transforms to direct matrix usage.
    // Pnm80():   NOT IMPLEMENTED: Prefer transforms to direct matrix usage. Outdated model.
    // P06e():    NOT IMPLEMENTED: Precession angles of uncertain relevance to amateur astronomy. Equinox.
    //            (Note: P06e() is out of order within SOFA's Library list-- it is where Po6e would be.)
    // Pom00():   NOT IMPLEMENTED: Prefer transforms to direct matrix usage.
    // Pr00():    NOT IMPLEMENTED: Precession rate of uncertain relevance to amateur astronomy.
    // Prec76():  NOT IMPLEMENTED: Precession model of low relevance to amateur astronomy. Outdated model.
    // S00():     NOT IMPLEMENTED: Outdated model.
    // S00a():    NOT IMPLEMENTED: Outdated model.
    // S00b():    NOT IMPLEMENTED: Outdated model.
    // S06():     NOT IMPLEMENTED: CIP X,Y coordinates typically not known.
    
    /// <summary> Return the CIO locator s, positioning the CIO on the equator of the
    /// Celestial Intermediate Pole, using the IAU 2006 precession and IAU 2000A nutation models.</summary>
    /// <param name="date1">Part 1 of Julian Date, TT</param>
    /// <param name="date2">Part 2 of Julian Date, TT</param>
    /// <returns>The CIO locator s, the difference between Right Ascension of ascending note of CIP equator
    /// in GCRS vs. (CIP,CIO) systems. Small fraction of 1 arcsecond in years 1900-2100 (radians)</returns>
    public static double S06a(double date1, double date2) {
        var s = S_S06a(date1, date2);
        return s; }
    [DllImport(DllFilename, EntryPoint = "iauS06a", CallingConvention = CallingConvention.Cdecl)]
    static extern double S_S06a(double date1, double date2);
    
    // Sp00():   NOT IMPLEMENTED: TIO locator of uncertain relevance to amateur astronomy.
    
    /// <summary> Sofa.Xy06(): Return X,Y coordinates of the Celestial Intermediate Pole,
    /// using a series based on IAU 2006 precession and IAU 2000A nutation models.</summary>
    /// <param name="date1">Part 1 of Julian Date, TT</param>
    /// <param name="date2">Part 2 of Julian Date, TT</param>
    /// <returns>2-tuple of doubles:  x: CIP X coordinate, y: CIP Y coordinate</returns>
    public static Tuple<double, double> Xy06(double date1, double date2) {
        double x = 0, y = 0;
        S_Xy06(date1, date2, ref x, ref y);
        return Tuple.Create(x, y); }
    [DllImport(DllFilename, EntryPoint = "iauXy06", CallingConvention = CallingConvention.Cdecl)]
    static extern void S_Xy06(double date1, double date2, ref double x, ref double y);
    
    // Xys00a(): Prefer Xys06a().
    // Xys00b(): Prefer Xys06a().
    
    /// <summary> Sofa.Xys06a():  For a given TT date, compute the X,Y coordinates of the
    /// Celestial Intermediate Pole and the CIO locator s, using the IAU 2006 precession
    /// and IAU 2000A nutation models.</summary>
    /// <param name="date1">Part 1 of Julian Date, TT</param>
    /// <param name="date2">Part 2 of Julian Date, TT</param>
    /// <returns> 3-Tuple of doubles: x: Celestial Intermediate Pole
    ///                               y: Celestial Intermediate Pole
    ///                               s: the CIO locator s</returns>
    public static Tuple<double, double, double> Xys06a(double date1, double date2) {
        double x = 0, y = 0, s = 0;
        S_Xys06a(date1, date2, ref x, ref y, ref s);
        return Tuple.Create(x, y, s); }
    [DllImport(DllFilename, EntryPoint = "iauXys06a", CallingConvention = CallingConvention.Cdecl)]
    static extern void S_Xys06a(double date1, double date2, ref double x, ref double y, ref double s);
    
    #endregion  //  Precession/Nutation/Polar Motion.


    #region "Astronomy: Star Catalog Conversions"
    /*######### Astronomy: Star Catalog Conversions ######################################################*/

    // All 9 functions in Star Catalog Conversions section:  NOT IMPLEMENTED.
    
    #endregion  //  Star Catalog Conversions.


    #region "Vector/Matrix: Initialization"
    /*######### Vector/Matrix: Initialization ############################################################*/
    // Zp():   NOT IMPLEMENTED. Implement in C# if needed.
    // Zr():   NOT IMPLEMENTED. Implement in C# if needed.
    // Ir():   NOT IMPLEMENTED. Implement in C# if needed.
    // Zpv():  NOT IMPLEMENTED. Implement in C# if needed.

    #endregion


    #region "Vector/Matrix: Copy/Extend/Extract"
    /*######### Vector/Matrix: Copy/Extend/Extract #######################################################*/
    // Cp():     NOT IMPLEMENTED. Implement in C# if needed.
    // Cr():     NOT IMPLEMENTED. Implement in C# if needed.
    // Cpv():    NOT IMPLEMENTED. Implement in C# if needed.
    // P2pv():   NOT IMPLEMENTED. Implement in C# if needed.
    // Pv2p():   NOT IMPLEMENTED. Implement in C# if needed.

    #endregion


    #region "Vector/Matrix: Build Rotations"
    /*######### Vector/Matrix: Build Rotations ###########################################################*/

    /// <summary> Sofa.Rx(): Rotate a rotation matrix about the x−axis.</summary>
    /// <param name="phi">Newly applied rotation angle (radians)</param>
    /// <param name="r">Rotation matrix before new rotation</param>
    /// <returns>same r: Rotation matrix after new rotation</returns>
    public static double[,] Rx(double phi, double[,] r) {
        S_Rx(phi, r);
        return r; }
    [DllImport(DllFilename, EntryPoint = "iauRx", CallingConvention = CallingConvention.Cdecl)]
    static extern void S_Rx(double phi, [In, Out] double[,] r);


    /// <summary> Sofa.Ry(): Rotate a rotation matrix about the y−axis.</summary>
    /// <param name="theta">Newly applied rotation angle (radians)</param>
    /// <param name="r">Rotation matrix before new rotation</param>
    /// <returns>same r: Rotation matrix after new rotation</returns>
    public static double[,] Ry(double theta, double[,] r) {
        S_Ry(theta, r);
        return r; }
    [DllImport(DllFilename, EntryPoint = "iauRy", CallingConvention = CallingConvention.Cdecl)]
    static extern void S_Ry(double theta, [In, Out] double[,] r);


    /// <summary> Sofa.Rz(): Rotate a rotation matrix about the z−axis.</summary>
    /// <param name="psi">Newly applied rotation angle (radians)</param>
    /// <param name="r">Rotation matrix before new rotation</param>
    /// <returns>same r: Rotation matrix after new rotation</returns>
    public static double[,] Rz(double psi, double[,] r) {
        S_Rz(psi, r);
        return r; }
    [DllImport(DllFilename, EntryPoint = "iauRz", CallingConvention = CallingConvention.Cdecl)]
    static extern void S_Rz(double psi, [In, Out] double[,] r);

    #endregion


    #region "Vector/Matrix: Spherical/Cartesian Conversions"
    /*######### Vector/Matrix: Spherical/Cartesian Conversions ###########################################*/
    
    /// <summary> Sofa.S2c(): Convert spherical coordinates to Cartesian unit vector
    /// (direction cosines).</summary>
    /// <param name="theta">Longitude angle (radians)</param>
    /// <param name="phi">Latitude angle (radians)</param>
    /// <returns>Unit vector result (direction cosines) as double[3]</returns>
    public static double[] S2c(double theta, double phi) {
        var c = new double[3] {0, 0, 0};
        S_S2c(theta, phi, c);
        return c; }
    [DllImport(DllFilename, EntryPoint = "iauS2c", CallingConvention = CallingConvention.Cdecl)]
    static extern void S_S2c(double theta, double phi, [In, Out] double[] c);
    
    
    /// <summary>Sofa.C2s(): Convert position vector to spherical coordinates.</summary>
    /// <param name="p">Position vector (double[3]) of any magnitude, only its direction is used.</param>
    /// <returns>2-Tuple:  theta (double): Longitude angle (radians)
    ///                    phi (double): Latitude angle (radians)</returns>
    public static Tuple<double, double> C2s(double[] p) {
        double theta = 0, phi = 0;
        S_C2s(p, ref theta, ref phi);
        return Tuple.Create(theta, phi); }
    [DllImport(DllFilename, EntryPoint = "iauC2s", CallingConvention = CallingConvention.Cdecl)]
    static extern double S_C2s(double[] p, ref double theta, ref double phi);


    /// <summary> Sofa.S2p(): Convert spherical polar coordinates to p−vector.</summary>
    /// <param name="theta">Longitude angle (radians)</param>
    /// <param name="phi">Latitude angle (radians)</param>
    /// <param name="r">Radial distance</param>
    /// <returns>Cartesian coordinates.</returns>
    public static double[] S2p(double theta, double phi, double r) {
        var p = new double[3] {0, 0, 0};
        S_S2p(theta, phi, r, p);
        return p; }
    [DllImport(DllFilename, EntryPoint = "iauS2p", CallingConvention = CallingConvention.Cdecl)]
    static extern void S_S2p(double theta, double phi, double r, [In, Out] double[] p);


    /// <summary>Sofa.P2s(): Convert position vector to spherical polar coordinates.</summary>
    /// <param name="p">Position vector (double[3]) of any magnitude</param>
    /// <returns>3-Tuple:  theta (double): Longitude angle (radians)
    ///                    phi (double): Latitude angle (radians)
    ///                    r (double): Radial distance</returns>
    public static Tuple<double, double, double> P2s(double[] p) {
        double theta = 0, phi = 0, r = 0;
        S_P2s(p, ref theta, ref phi, ref r);
        return Tuple.Create(theta, phi, r); }
    [DllImport(DllFilename, EntryPoint = "iauP2s", CallingConvention = CallingConvention.Cdecl)]
    static extern void S_P2s(double[] p, ref double theta, ref double phi, ref double r);


    /// <summary> Sofa.S2pv(): Convert position/velocity from spherical to Cartesian coordinates.
    /// The inverse of function Sofa.Pv2s().</summary>
    /// <param name="theta">Longitude angle (radians)</param>
    /// <param name="phi">Latitude angle (radians)</param>
    /// <param name="r">Radial distance</param>
    /// <param name="td">Rate of change of theta</param>
    /// <param name="pd">Rate of change of phi</param>
    /// <param name="rd">Rate of change of radius</param>
    /// <returns>Position/velocity vector as double[2,3]</returns>
    public static double[,] S2pv(double theta, double phi, double r, double td, double pd, double rd) {
        var pv = new double[2, 3] {{0, 0, 0}, {0, 0, 0}};
        S_S2pv(theta, phi, r, td, pd, rd, pv);
        return pv; }
    [DllImport(DllFilename, EntryPoint = "iauS2pv", CallingConvention = CallingConvention.Cdecl)]
    static extern void S_S2pv(double theta, double phi, double r,
        double td, double pd, double rd, [In, Out] double[,] pv);


    /// <summary>Sofa.Pv2s(): Convert position/velocity from Cartesian to spherical coordinates
    /// The inverse of function Sofa.S2pv().</summary>
    /// <param name="pv">Position/velocity vectors as double[2,3]</param>
    /// <returns>6-Tuple:  theta (double): Longitude angle (radians)
    ///                    phi (double): Latitude angle (radians)
    ///                    r (double): Radial distance
    ///                    td (double): Rate of change of theta
    ///                    pd (double): Rate of change of phi
    ///                    rd (double): Rate of change of radius</returns>
    public static Tuple<double, double, double, double, double, double> Pv2s(double[,] pv) {
        double theta = 0, phi = 0, r = 0, td = 0, pd = 0, rd = 0;
        S_Pv2s(pv, ref theta, ref phi, ref r, ref td, ref pd, ref rd);
        return Tuple.Create(theta, phi, r, td, pd, rd); }
    [DllImport(DllFilename, EntryPoint = "iauPv2s", CallingConvention = CallingConvention.Cdecl)]
    static extern void S_Pv2s(double[,] pv, ref double theta, ref double phi, ref double r,
        ref double td, ref double pd, ref double rd);

    #endregion //  Vector/Matrix: Spherical/Cartesian Conversions


    #region "Vector/Matrix: Operations on Vectors"
    /*######### Vector/Matrix: Operations on Vectors #####################################################*/

    /*##### For position (p-) vectors and rotation (r-) matrices: #####*/
    // Ppp():     NOT IMPLEMENTED. Implement in C# if needed.
    // Pmp():     NOT IMPLEMENTED. Implement in C# if needed.
    // Ppsp():    NOT IMPLEMENTED. Implement in C# if needed.
    // Pdp():     NOT IMPLEMENTED. Implement in C# if needed.
    // Pxp():     NOT IMPLEMENTED. Implement in C# if needed.
    // Pm():      NOT IMPLEMENTED. Implement in C# if needed.
    // Pn():      NOT IMPLEMENTED. Implement in C# if needed.
    // Sxp():     NOT IMPLEMENTED. Implement in C# if needed.

    /*##### For position-velocity (pv-) vectors: #####*/
    // Pvppv():   NOT IMPLEMENTED. Implement in C# if needed.
    // Pvmpv():   NOT IMPLEMENTED. Implement in C# if needed.
    // Pvdpv():   NOT IMPLEMENTED. Implement in C# if needed.
    // Pvxpv():   NOT IMPLEMENTED. Implement in C# if needed.
    // Pvm():     NOT IMPLEMENTED. Implement in C# if needed.
    // Sxpv():    NOT IMPLEMENTED. Implement in C# if needed.
    // S2xpv():   NOT IMPLEMENTED. Implement in C# if needed.
    // Pvu():     NOT IMPLEMENTED. Implement in C# if needed.
    // Pvup():    NOT IMPLEMENTED. Implement in C# if needed.

    #endregion //  Vector/Matrix: Operations on Vectors


    #region "Vector/Matrix: Operations on Matrices"
    /*######### Vector/Matrix: Operations on Matrices ####################################################*/

    /// <summary> Sofa.Rxr(): Multiply two rotation matrices.
    /// Same array may be used for any of the arguments.
    /// NB: not commutative, Rxr(a, b) != Rxr(b, a).</summary>
    /// <param name="a">First rotation matrix</param>
    /// <param name="b">Second rotation matrix</param>
    /// <returns>New rotation matrix [3,3]: a * b</returns>
    public static double[,] Rxr(double[,] a, double[,] b) {
        var atb = new double[3, 3] {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
        S_Rxr(a, b, atb);
        return atb; }
    [DllImport(DllFilename, EntryPoint = "iauRxr", CallingConvention = CallingConvention.Cdecl)]
    static extern void S_Rxr(double[,] a, double[,] b, [In, Out] double[,] atb);


    /// <summary> Sofa.Tr(): Return transpose of a rotation matrix.</summary>
    /// <param name="r">Rotation matrix before transposition</param>
    /// <returns>Transposed rotation matrix</returns>
    public static double[,] Tr(double[,] r) {
        var rt = new double[3, 3] {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
        S_Tr(r, rt);
        return rt; }
    [DllImport(DllFilename, EntryPoint = "iauTr", CallingConvention = CallingConvention.Cdecl)]
    static extern void S_Tr(double[,] r, [In, Out] double[,] rt);

    #endregion  //  Operations on Matrices.


    #region "Vector/Matrix: Matrix-Vector Products"
    /*######### Vector/Matrix: Matrix-Vector Products ####################################################*/
    
    /// <summary>Sofa.Rxp(): Multiply a position vector by a rotation matrix,
    /// return the rotated vector.</summary>
    /// <param name="r">Rotation matrix</param>
    /// <param name="p">Position vector before rotation</param>
    /// <returns>Position vector after rotation, i.e., r * p</returns>
    public static double[] Rxp(double[,] r, double[] p) {
        var rp = new double[3] {0, 0, 0};
        S_Rxp(r, p, rp);
        return rp; }
    [DllImport(DllFilename, EntryPoint = "iauRxp", CallingConvention = CallingConvention.Cdecl)]
    static extern void S_Rxp(double[,] r, double[] p, [In, Out] double[] rp);

    // Trxp():     NOT IMPLEMENTED. Transpose form; implement later if needed. 

    /// <summary> Sofa.Rxpv(): Multiply a position/velocity (pv-) vector by a rotation (r-) matrix,
    /// returning the rotated pv-vector.</summary>
    /// <param name="r">Rotation (r-) matrix</param>
    /// <param name="pv">Position/velocity (pv-) vector</param>
    /// <returns>Rotated position/velocity (pv-) vector</returns>
    public static double[,] Rxpv(double[,] r, double[,] pv) {
        var rpv = new double[2,3] {{0, 0, 0}, {0, 0, 0}};
        S_Rxpv(r, pv, rpv);
        return rpv; }
    [DllImport(DllFilename, EntryPoint = "iauRxpv", CallingConvention = CallingConvention.Cdecl)]
    static extern void S_Rxpv(double[,] r, double[,] pv, [In, Out] double[,] rpv);

    // Trxpv():    NOT IMPLEMENTED. Transpose form; implement later if needed. 

    #endregion //  Matrix-Vector Products


    #region "Vector/Matrix: Separation and Position Angle"
    /*######### Vector/Matrix: Separation and Position Angle #############################################*/

    /// <summary>Sofa.Sepp(): Angular separation between two p−vectors.</summary>
    /// <param name="a">First p−vector (not necessarily unit length)</param>
    /// <param name="b">Second p−vector (not necessarily unit length)</param>
    /// <returns>Angular separation (radians, always in range [0, pi])</returns>
    public static double Sepp(double[] a, double[] b) {
        var separation = S_Sepp(a, b);
        return separation; }
    [DllImport(DllFilename, EntryPoint = "iauSepp", CallingConvention = CallingConvention.Cdecl)]
    static extern double S_Sepp(double[] a, double[] b);


    /// <summary>Sofa.Seps(): Angular separation between two sets of spherical coordinates.</summary>
    /// <param name="al">First longitude/RA (radians)</param>
    /// <param name="ap">First latitude/Declination (radians)</param>
    /// <param name="bl">Second longitude/RA (radians)</param>
    /// <param name="bp">Second latitude/Declination (radians)</param>
    /// <returns>Angular separation (radians)</returns>
    public static double Seps(double al, double ap, double bl, double bp) {
        var separation = S_Seps(al, ap, bl, bp);
        return separation; }
    [DllImport(DllFilename, EntryPoint = "iauSeps", CallingConvention = CallingConvention.Cdecl)]
    static extern double S_Seps(double al, double ap, double bl, double bp);

    // Pap():    NOT IMPLEMENTED. Implement later if position angles needed.
    // Pas():    NOT IMPLEMENTED. Implement later if position angles needed.

    #endregion //  Separation and Position Angle


    #region "Vector/Matrix: Rotation Vectors"
    /*######### Vector/Matrix: Rotation Vectors ##########################################################*/

    // Rm2v():   NOT IMPLEMENTED. Implement later if needed.
    // Rv2m():   NOT IMPLEMENTED. Implement later if needed.

    #endregion  //  Rotation Vectors.


    #region "Vector/Matrix: Operations on Angles"
    /*######### Vector/Matrix: Operations on Angles ######################################################*/

    /// <summary> Sofa.A2af(): Decompose radians into degrees, arcminutes, arcseconds, fraction.</summary>
    /// <param name="ndp">Resolution, as negative 10-exponent of last decimal place,
    /// e.g., -2 for rounding to hundreds, +1 for rounding to 1/10, +6 for rounding to millionths</param>
    /// <param name="angle">Input angle (radians). May exceed 2*pi.
    /// NB: user is responsible for handling when near 360 degrees</param>
    /// <returns>2-Tuple:
    ///     sign (string of length one): '+' or '-'.
    ///     idmsf (int array of length 4):
    ///         degrees, arcminutes, arcseconds, fraction in the requested resolution.</returns>
    public static Tuple<string, int[]> A2af(int ndp, double angle) {
        char[] signChar = {'X'};
        var idmsf = new int[4] {0, 0, 0, 0};
        S_A2af(ndp, angle, signChar, idmsf);
        string signString = new string(signChar);
        return Tuple.Create(signString, idmsf); }
    [DllImport(DllFilename, EntryPoint = "iauA2af", CallingConvention = CallingConvention.Cdecl)]
    static extern void S_A2af(int ndp, double angle, [In, Out] char[] signChar, [In, Out] int[] idmsf);

    
    /// <summary> Sofa.A2tf(): Decompose radians into hours, minutes, seconds, fraction.</summary>
    /// <param name="ndp">Resolution, as negative 10-exponent of last decimal place,
    /// e.g., -2 for rounding to hundreds, +1 for rounding to 1/10, +6 for rounding to millionths</param>
    /// <param name="angle">Input angle (radians). May exceed 2*pi.
    /// NB: user is responsible for handling when near 360 degrees.</param>
    /// <returns>2-Tuple:
    ///     sign (string of length one): '+' or '-'.
    ///     ihmsf (int array of length 4):
    ///         hours, minutes, seconds, fraction in the requested resolution.</returns>
    public static Tuple<string, int[]> A2tf(int ndp, double angle) {
        char[] signChar = {'X'};
        var ihmsf = new int[4] {0, 0, 0, 0};
        S_A2tf(ndp, angle, signChar, ihmsf);
        string signString = new string(signChar);
        return Tuple.Create(signString, ihmsf); }
    [DllImport(DllFilename, EntryPoint = "iauA2tf", CallingConvention = CallingConvention.Cdecl)]
    static extern void S_A2tf(int ndp, double angle, [In, Out] char[] signChar, [In, Out] int[] ihmsf);

    
    /// <summary> Sofa.Af2a(): Convert degrees, arcminutes, arcseconds to radians.</summary>
    /// <param name="s">Sign [string], '-' means negative, any other means positive.</param>
    /// <param name="ideg">Degrees [int]</param>
    /// <param name="iamin">Arcminutes [int]</param>
    /// <param name="asec">Arcseconds [double]</param>
    /// <returns> 2-Tuple: status [int]: 0 -> OK
    ///                                  1 -> ideg outside [0-359] (degrees)
    ///                                  2 -> iamin outside [0-59] (minutes)
    ///                                  3 -> asec outside [0, 60) (arcseconds)
    ///                    rad: Angle (radians)</returns>
    public static Tuple<int, double> Af2a(string s, int ideg, int iamin, double asec) {
        double rad = 0;
        char sChar = s[0];
        var status = S_Af2a(sChar, ideg, iamin, asec, ref rad);
        return Tuple.Create(status, rad); }
    [DllImport(DllFilename, EntryPoint = "iauAf2a", CallingConvention = CallingConvention.Cdecl)]
    static extern int S_Af2a(char sChar, int ideg, int iamin, double asec, ref double rad);

    // Anp():   NOT IMPLEMENTED. Prefer AstroLib.Core.AstroMath.Wrap().
    // Anpm():  NOT IMPLEMENTED. Prefer AstroLib.Core.AstroMath.Wrap().

    /// <summary> Sofa. D2tf(): Decompose days to hours, minutes, seconds, fraction.</summary>
    /// <param name="ndp">Resolution, as negative 10-exponent of last decimal place,
    /// e.g., -2 for rounding to hundreds, +1 for rounding to 1/10, +6 for rounding to millionths.</param>
    /// <param name="days">Interval in days.</param>
    /// <returns>2-Tuple:
    ///     sign (string of length one): '+' or '-'.
    ///     ihmsf (int array of length 4):
    ///         hours, minutes, seconds, fraction in the requested resolution.</returns>
    public static Tuple<string, int[]> D2tf(int ndp, double days) {
        char[] signChar = {'X'};
        var ihmsf = new int[4] {0, 0, 0, 0};
        S_D2tf(ndp, days, signChar, ihmsf);
        string signString = new string(signChar);
        return Tuple.Create(signString, ihmsf); }
    [DllImport(DllFilename, EntryPoint = "iauD2tf", CallingConvention = CallingConvention.Cdecl)]
    static extern void S_D2tf(int ndp, double days, [In, Out] char[] signChar, [In, Out] int[] ihmsf);
    

    /// <summary> Sofa.Tf2a(): Convert hours, minutes, seconds to angle in radians.</summary>
    /// <param name="s">Sign [string], '-' means negative, any other means positive.</param>
    /// <param name="ihour">Hours (int)</param>
    /// <param name="imin">Minutes (int)</param>
    /// <param name="sec">Seconds (double)</param>
    /// <returns> 2-Tuple: status [int]: 0 -> OK
    ///                                  1 -> ideg outside [0-359] (degrees)
    ///                                  2 -> iamin outside [0-59] (minutes)
    ///                                  3 -> asec outside [0, 60) (arcseconds)
    ///                    rad: Angle (radians)</returns>
    public static Tuple<int, double> Tf2a(string s, int ihour, int imin, double sec) {
        double rad = 0;
        char ch = s[0];
        var status = S_Tf2a(ch, ihour, imin, sec, ref rad);
        return Tuple.Create(status, rad); }
    [DllImport(DllFilename, EntryPoint = "iauTf2a", CallingConvention = CallingConvention.Cdecl)]
    static extern int S_Tf2a(char ch, int ihour, int imin, double sec, ref double rad);


    /// <summary> Sofa.Tf2d(): Convert hours, minutes, seconds to days. </summary>
    /// <param name="s">Sign string, length 1], '-' means negative, any other means positive.</param>
    /// <param name="ihour">Hours (int)</param>
    /// <param name="imin">Minutes (int)</param>
    /// <param name="sec">Seconds (double)</param>
    /// <returns>Days (double)</returns>
    public static Tuple<int, double> Tf2d(string s, int ihour, int imin, double sec) {
        char ch = s[0];
        double days = 0;
        var status = S_Tf2d(ch, ihour, imin, sec, ref days);
        return Tuple.Create(status, days); }
    [DllImport(DllFilename, EntryPoint = "iauTf2d", CallingConvention = CallingConvention.Cdecl)]
    static extern int S_Tf2d(char ch, int ihour, int imin, double sec, ref double days);

    #endregion  //  Operations on Angles.
    
}  //  class Sofa.
