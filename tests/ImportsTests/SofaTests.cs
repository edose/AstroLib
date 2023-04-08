using System.Diagnostics.CodeAnalysis;
using AstroLib.Imports.Sofa;
using NUnit.Framework;
// ReSharper disable CompareOfFloatsByEqualityOperator
// ReSharper disable CommentTypo
// ReSharper disable InconsistentNaming
// ReSharper disable IdentifierTypo
// ReSharper disable InlineTemporaryVariable
// ReSharper disable JoinDeclarationAndInitializer

namespace AstroLib.ImportsTests2;

[TestFixture]
[SuppressMessage("Assertion", "NUnit2045:Use Assert.Multiple")]
public class SofaTests {
    private static double DegreesAsRadians(double deg) {
        return deg * (Math.PI / 180.0);
    }

    private static double HoursAsRadians(double hours) {
        return (hours * 15.0) * (Math.PI / 180.0);
    }

    private static double RadiansAsDegrees(double radians) {
        return radians * (180.0 / Math.PI);
    }

    private const double ra_yy_gem = 7.5770553; // hours
    private const double dec_yy_gem = 31.86944; // degrees
    private const double nms_longitude = -105.0 - (32.0 / 60.0) - (44.0 / 3600.0); // degrees
    private const double nms_latitude = 32.0 + (54.0 / 60.0) + (11.0 / 3600.0); // degrees
    private const double nms_height = 2228.0; //meters above WGS84
    private const double jd1 = 2400000.5, jd2 = 55000.0;
    private const double jdJ2000A = 2400000.5;
    private const double jdJ2000B = 2451545.0 - jdJ2000A;
    
    // values on date of UTC 2023-02-24:
    private const double ut1_utc = -0.0144; // UT1-UTC, "dut1".
    private const double tai_utc = 37.0;    // TAI-UTC, cumulative leap seconds, "dat".
    private const double tt_tai = 32.184;   // TT-TAI, by definition.
    // Common aliases:
    private const double dut1_ = ut1_utc;                  // UT1-UTC, classical "delta UT1".
    private const double dta_ = ut1_utc - tai_utc;         // UT1-TAI, "dta".
    private const double dat_ = tai_utc;                   // TAI-UTC, classical delta(AT)". 
    private const double dt_ = tt_tai + tai_utc - ut1_utc; // TT-UT1, classical "delta T".

    // [SetUp]
    // public void Setup() {}

#region "Astronomy/Calendars"
    
    [Test]
    public void Cal2jd_Tests() {
        // Cal2jd(): Gregorian calendar year, month day -> JD UTC.
        // Tests pass: 2023-02-25.
        var (djm0, djm) = Sofa.Cal2jd(2023, 2, 24);
        Assert.That(djm0, Is.TypeOf<double>());
        Assert.That(djm, Is.TypeOf<double>());
        Assert.That(djm0, Is.EqualTo(2400000.5));
        Assert.That(djm, Is.EqualTo(59999));
    }

    [Test]
    public void Jd2cal_Tests() {
        // JD2cal(): JD UTC -> Gregorian calendar year, month day.
        // Tests pass: 2023-02-25.
        var (iy, im, id, fd) = Sofa.Jd2cal(2450000.5, 9999.125);
        Assert.That(iy, Is.TypeOf<int>());
        Assert.That(im, Is.TypeOf<int>());
        Assert.That(id, Is.TypeOf<int>());
        Assert.That(fd, Is.TypeOf<double>());
        Assert.That((iy, im, id), Is.EqualTo((2023, 2, 24)));
        Assert.That(fd, Is.EqualTo(0.125));
    }

    #endregion

#region "Astronomy/Astrometry"

    [Test]
    public void Atcc13_Tests() {
        // Atcc13(): ICRS J2000 -> ICRS astrometric at given date.
        // Tests pass: 2023-02-25.
        var rc = HoursAsRadians(6.0);
        var dc = DegreesAsRadians(35.0);
        var (ra, da) = Sofa.Atcc13(
            rc, dc,
            0.000_001, -0.000_002, 0, 0,
            2400000.5, 55000);
        Assert.That(ra, Is.TypeOf<double>());
        Assert.That(da, Is.TypeOf<double>());
        Assert.That(ra, Is.EqualTo(1.570_805_787).Within(0.000_000_001));
        Assert.That(da, Is.EqualTo(0.610_846_317).Within(0.000_000_001));
    }

    [Test]
    public void Atci13_Tests() {
        // Atci13(): ICRS J2000.0 -> CIRS at given date.
        // Tests pass: 2023-02-25.
        var rc = HoursAsRadians(ra_yy_gem);
        var dc = DegreesAsRadians(dec_yy_gem);
        var (ri, di, eo) = Sofa.Atci13(
            rc, dc, 0.000_000, -0.000_000, 0, 0, jd1, jd2);
        Assert.That(ri, Is.TypeOf<double>());
        Assert.That(di, Is.TypeOf<double>());
        Assert.That(eo, Is.TypeOf<double>());
        var ri_expected = DegreesAsRadians(113.681_104_28);
        var di_expected = DegreesAsRadians(31.849_973_63);
        Assert.That(ri, Is.EqualTo(ri_expected).Within(0.000_001));
        Assert.That(di, Is.EqualTo(di_expected).Within(0.000_001));
        Assert.That(eo, Is.EqualTo(-0.002_177_743).Within(0.000_000_001));
    }

    [Test]
    public void Atco13_Tests() {
        // ICRS RA,Dec to observables at given site and time.
        // Tests pass: 2023-02-25 (and match TheSkyX results).
        var rc = HoursAsRadians(7.5770553); // YY Gem, ICRS
        var dc = DegreesAsRadians(31.86944);
        const double pr = 0.0;
        const double pd = 0.0;
        const double px = 0.0;
        const double rv = 0.0;
        const double utc1 = 2400000.5;
        const double utc2 = 55000;
        const double dut1 = -0.0099;
        var elong = DegreesAsRadians(nms_longitude); // NMS
        var phi = DegreesAsRadians(nms_latitude);
        const double height = 2228;
        const double xp = 0;
        const double yp = 0;
        const double phpa = 800;
        const double tc = 0;
        const double rh = 0.20;
        const double wl = 0.60;
        var (aob, zob, hob, dob, rob, eo) =
            Sofa.Atco13(rc, dc, pr, pd, px, rv, utc1, utc2,
                dut1, elong, phi, height, xp, yp, phpa, tc, rh, wl);
        Assert.That(aob, Is.TypeOf<double>());
        Assert.That(zob, Is.TypeOf<double>());
        Assert.That(hob, Is.TypeOf<double>());
        Assert.That(dob, Is.TypeOf<double>());
        Assert.That(rob, Is.TypeOf<double>());
        Assert.That(eo, Is.TypeOf<double>());
        Assert.That(aob, Is.EqualTo(DegreesAsRadians(281.680)).Within(0.001));
        Assert.That(zob, Is.EqualTo(DegreesAsRadians(90.0 - 50.624)).Within(0.001));
        Assert.That(hob, Is.EqualTo(DegreesAsRadians(47.005)).Within(0.001));
        Assert.That(dob, Is.EqualTo(DegreesAsRadians(31.853)).Within(0.001));
        Assert.That(rob, Is.EqualTo(DegreesAsRadians(113.694)).Within(0.001));
        Assert.That(eo, Is.EqualTo(-0.002_178).Within(0.000_001));
    }

    [Test]
    public void Atic13_Tests() {
        // Atic13(): CIRS RA,Dec (geocentric) -> ICRS astrometric.
        // Tests pass: 2023-02-25 (and match astropy results).
        var ri = HoursAsRadians(6.0);
        var di = DegreesAsRadians(35.0);
        const double date1 = 2400000.5;
        const double date2 = 55000.0;
        var (rc, dc, eo) = Sofa.Atic13(ri, di, date1, date2);
        Assert.That(rc, Is.TypeOf<double>());
        Assert.That(dc, Is.TypeOf<double>());
        Assert.That(eo, Is.TypeOf<double>());
        Assert.That(rc, Is.EqualTo(DegreesAsRadians(89.968_864_04)).Within(0.000_001));
        Assert.That(dc, Is.EqualTo(DegreesAsRadians(34.998_734_45)).Within(0.000_001));
        Assert.That(eo, Is.EqualTo(-0.002_177_743).Within(0.000_000_001));
    }

    [Test]
    public void Atio13_Tests() {
        // Atic13(): CIRS RA,Dec (geocentric) -> observables at given site and time.
        // Tests pass: 2023-02-25.
        var rc = HoursAsRadians(7.5770553); // YY Gem, ICRS
        var dc = DegreesAsRadians(31.86944);
        const double pr = 0.0;
        const double pd = 0.0;
        const double px = 0.0;
        const double rv = 0.0;
        const double dut1 = -0.0099;
        var elong = DegreesAsRadians(-105.0 - (32.0 / 60.0) - (44.0 / 3600.0)); // NMS
        var phi = DegreesAsRadians(32.0 + (54.0 / 60.0) + (11.0 / 3600.0));
        const double height = 2228;
        const double xp = 0;
        const double yp = 0;
        const double phpa = 800;
        const double tc = 0;
        const double rh = 0.20;
        const double wl = 0.60;
        // RA,Dec ICRS -> CIRS via Atci13() which is tested above:
        var (ri, di, eo) = Sofa.Atci13(
            rc, dc, pr, pd, px, rv, jd1, jd2);
        // Get Observables from CIRS via Atio13():
        var (aob, zob, hob, dob, rob) = Sofa.Atio13(ri, di,
            jd1, jd2, dut1, elong, phi, height, xp, yp, phpa, tc, rh, wl);
        Assert.That(aob, Is.EqualTo(DegreesAsRadians(281.680)).Within(0.001));
        Assert.That(zob, Is.EqualTo(DegreesAsRadians(90.0 - 50.624)).Within(0.001));
        Assert.That(hob, Is.EqualTo(DegreesAsRadians(47.005)).Within(0.001));
        Assert.That(dob, Is.EqualTo(DegreesAsRadians(31.853)).Within(0.001));
        Assert.That(rob, Is.EqualTo(DegreesAsRadians(113.694)).Within(0.001));
        Assert.That(eo, Is.EqualTo(-0.002_178).Within(0.000_001));
    }

    [Test]
    public void Atoc13_Tests() {
        // Atoc13(): observables at given site and time -> ICRS astrometric.
        // Tests pass: 2023-02-25.
        const string type = "A";
        var ob1 = DegreesAsRadians(281.680);
        var ob2 = DegreesAsRadians(90.0 - 50.624);
        const double dut1 = -0.0099;
        var elong = DegreesAsRadians(nms_longitude); // NMS
        var phi = DegreesAsRadians(nms_latitude);
        const double height = 2228;
        const double xp = 0;
        const double yp = 0;
        const double phpa = 800;
        const double tc = 0;
        const double rh = 0.20;
        const double wl = 0.60;
        var (rc, dc) = Sofa.Atoc13(type, ob1, ob2, jd1, jd2,
            dut1, elong, phi, height, xp, yp, phpa, tc, rh, wl);
        var rc_expected = HoursAsRadians(7.5770553); // YY Gem, ICRS
        var dc_expected = DegreesAsRadians(31.86944);
        Assert.That(rc, Is.EqualTo(rc_expected).Within(0.001));
        Assert.That(dc, Is.EqualTo(dc_expected).Within(0.001));
    }

    [Test]
    public void Atoi13_Tests() {
        // Atoi13(): observables at given site and time -> CIRS.
        // Tests pass: 2023-02-25.
        const string type = "A";
        var ob1 = DegreesAsRadians(281.680);
        var ob2 = DegreesAsRadians(90.0 - 50.624);
        const double dut1 = -0.0099;
        var elong = DegreesAsRadians(-105.0 - (32.0 / 60.0) - (44.0 / 3600.0)); // NMS
        var phi = DegreesAsRadians(32.0 + (54.0 / 60.0) + (11.0 / 3600.0));
        const double height = 2228;
        const double xp = 0;
        const double yp = 0;
        const double phpa = 800;
        const double tc = 0;
        const double rh = 0.20;
        const double wl = 0.60;
        var (ri, di) = Sofa.Atoi13(type, ob1, ob2, jd1, jd2,
            dut1, elong, phi, height, xp, yp, phpa, tc, rh, wl);
        // Expected values via Atoc13() and Atci13(), both tested above:
        var (rc, dc) = Sofa.Atoc13(type, ob1, ob2, jd1, jd2,
            dut1, elong, phi, height, xp, yp, phpa, tc, rh, wl);
        const double pr = 0;
        const double pd = 0;
        const double px = 0;
        const double rv = 0;
        var (riExpected, diExpected, eo) = Sofa.Atci13(
            rc, dc, pr, pd, px, rv, jd1, jd2);
        Assert.That(ri, Is.EqualTo(riExpected).Within(0.001));
        Assert.That(di, Is.EqualTo(diExpected).Within(0.001));
    }

    [Test]
    public void Pmsafe_Tests() {
        // Pmsafe(): Update star catalog data for proper motion etc.
        // Tests pass: 2023-02-25.
        const double pmr1 = 0.000_001; // radians/year
        const double pmd1 = -0.000_002;
        const double px1 = 0.0;
        const double rv1 = 0.0;
        const double jdBeforeA = jdJ2000A;
        const double jdBeforeB = jdJ2000B;
        const double yearsElapsed = 4.5;
        const double jdAfterA = jdBeforeA;
        const double jdAfterB = jdBeforeB + yearsElapsed * 365.25;
        var raBeforeAsRadians = HoursAsRadians(ra_yy_gem);
        var decBeforeAsRadians = DegreesAsRadians(dec_yy_gem);
        var (ra2, dec2, pmr2, pmd2, px2, rv2) =
            Sofa.Pmsafe(raBeforeAsRadians, decBeforeAsRadians,
                pmr1, pmd1, px1, rv1, jdBeforeA, jdBeforeB, jdAfterA, jdAfterB);
        Assert.That(raBeforeAsRadians, Is.EqualTo(1.983_668_439).Within(0.000_000_100)); // input check
        Assert.That(decBeforeAsRadians, Is.EqualTo(0.556_226_659).Within(0.000_000_100)); // "
        Assert.That(ra2, Is.EqualTo(1.983_672_939).Within(0.000_000_100)); // result
        Assert.That(dec2, Is.EqualTo(0.556_217_659).Within(0.000_000_100)); // result
    }

    [Test]
    public void Pvtob_Tests() {
        // Pvtob(): Position and velocity (CIRS) of a terrestrial observing station.
        // Tests pass: 2023-02-25.
        const double xp = 0.0, yp = 0.0, sp = 0.0; // approximation
        double theta = 0.0; // Earth rotation angle
        var pv = new double[2, 3];
        pv = Sofa.Pvtob(DegreesAsRadians(nms_longitude), DegreesAsRadians(nms_latitude),
            nms_height, xp, yp, sp, theta);
        Assert.That(pv[0, 0], Is.EqualTo(-1437_091.805).Within(0.001)); // meters
        Assert.That(pv[0, 1], Is.EqualTo(-5166031.908).Within(0.001));
        Assert.That(pv[0, 2], Is.EqualTo(3446147.028).Within(0.001));
        Assert.That(pv[1, 0], Is.EqualTo(376.7130).Within(0.001)); // m/s
        Assert.That(pv[1, 1], Is.EqualTo(-104.7944).Within(0.001));
        Assert.That(pv[1, 2], Is.EqualTo(0.0000).Within(0.001));
        Console.WriteLine();
        pv = Sofa.Pvtob(DegreesAsRadians(nms_longitude), DegreesAsRadians(nms_latitude),
            nms_height, xp, yp, sp, 0.5 * Math.PI);
        Assert.That(pv[1, 0], Is.EqualTo(104.7944).Within(0.001)); // m/s
        Assert.That(pv[1, 1], Is.EqualTo(376.7130).Within(0.001));
        Assert.That(pv[1, 2], Is.EqualTo(0.0000).Within(0.001));
    }

    [Test]
    public void Refco_Tests() {
        // Refco(): Determine the constants in the atmospheric refraction model.
        // Tests pass: 2023-02-25.
        const double phpa = 1005; // example from SOFA manual.
        const double tc = 7;
        const double rh = 0.80;
        const double wl = 0.5740;
        var (refa, refb) = Sofa.Refco(phpa, tc, rh, wl);
        Assert.That(refa, Is.EqualTo(0.000_282_371).Within(0.000_000_1));
        Assert.That(refb, Is.EqualTo(-0.000_000_312).Within(0.000_000_001));
        var Z = DegreesAsRadians(60.0); // zenith distance
        var tanZ = Math.Tan(Z);
        var dZ = refa * tanZ + refb * Math.Pow(tanZ, 3.0);
        var dZ_arcseconds = 3600.0 * RadiansAsDegrees(dZ);
        Assert.That(dZ_arcseconds, Is.EqualTo(100.546).Within(0.01)); // example from SOFA manual.
    }

    #endregion

#region "Astronomy/Ephemerides"

    [Test]
    public void Epv00_Tests() {
        // Epv00(): Earth position and velocity, heliocentric and barycentric, BCRS.
        // Tests pass: 2023-02-25.
        var pvh = new double[2, 3];
        var pvb = new double[2, 3];
        double date1, date2;
        (date1, date2) = Sofa.Cal2jd(2023, 2, 24);
        (pvh, pvb) = Sofa.Epv00(date1, date2);
        Assert.That(pvh[0, 0], Is.EqualTo(-0.895_206).Within(0.000_001));
        Assert.That(pvh[1, 1], Is.EqualTo(-0.014_347).Within(0.000_001));
        Assert.That(pvb[0, 0], Is.EqualTo(-0.904_191).Within(0.000_001));
        Assert.That(pvb[1, 2], Is.EqualTo(-0.006_223).Within(0.000_001));

        (date1, date2) = Sofa.Cal2jd(2023, 5, 24);
        (pvh, pvb) = Sofa.Epv00(date1, date2);
        Assert.That(pvh[0, 0], Is.EqualTo(-0.470_967).Within(0.000_001));
        Assert.That(pvh[1, 1], Is.EqualTo(-0.007_401).Within(0.000_001));
        Assert.That(pvb[0, 0], Is.EqualTo(-0.479_765).Within(0.000_001));
        Assert.That(pvb[1, 2], Is.EqualTo(-0.003_212).Within(0.000_001));
    }

    [Test]
    public void Moon98_Tests() {
        // Moon98(): Moon position and velocity, geocentric, GCRS.
        // Tests pass: 2023-02-26.
        var pv = new double[2, 3];
        double date1, date2;
        (date1, date2) = Sofa.Cal2jd(2023, 2, 24);
        pv = Sofa.Moon98(date1, date2);
        Assert.That(pv[0, 0], Is.EqualTo(0.002_278_0).Within(0.000_000_1));
        Assert.That(pv[0, 1], Is.EqualTo(0.000_984_6).Within(0.000_000_1));
        Assert.That(pv[0, 2], Is.EqualTo(0.000_381_7).Within(0.000_000_1));
        Assert.That(pv[1, 0], Is.EqualTo(-0.000_213_8).Within(0.000_000_1));
        Assert.That(pv[1, 1], Is.EqualTo(0.000_495_6).Within(0.000_000_1));
        Assert.That(pv[1, 2], Is.EqualTo(0.000_271_3).Within(0.000_000_1));

        (date1, date2) = Sofa.Cal2jd(2023, 3, 3);
        pv = Sofa.Moon98(date1, date2);
        Assert.That(pv[0, 0], Is.EqualTo(-0.001_017_8).Within(0.000_000_1));
        Assert.That(pv[0, 1], Is.EqualTo(0.002_202_5).Within(0.000_000_1));
        Assert.That(pv[0, 2], Is.EqualTo(0.001_210_1).Within(0.000_000_1));
        Assert.That(pv[1, 0], Is.EqualTo(-0.000_518_9).Within(0.000_000_1));
        Assert.That(pv[1, 1], Is.EqualTo(-0.000_194_3).Within(0.000_000_1));
        Assert.That(pv[1, 2], Is.EqualTo(-0.000_071_5).Within(0.000_000_1));
    }

    [Test]
    public void Plan94_Tests() {
        // Plan94(): Planet position and velocity, geocentric, GCRS.
        // Tests pass: 2023-02-26.
        var pv = new double[2, 3];
        double date1, date2;
        (date1, date2) = Sofa.Cal2jd(2023, 2, 24);
        pv = Sofa.Plan94(date1, date2, 1); // Mercury
        Assert.That(pv[0, 0], Is.EqualTo(0.079_841_5).Within(0.000_000_1));
        Assert.That(pv[0, 1], Is.EqualTo(-0.393_789_0).Within(0.000_000_1));
        Assert.That(pv[0, 2], Is.EqualTo(-0.218_637_2).Within(0.000_000_1));
        Assert.That(pv[1, 0], Is.EqualTo(0.022_055_8).Within(0.000_000_1));
        Assert.That(pv[1, 1], Is.EqualTo(0.006_446_1).Within(0.000_000_1));
        Assert.That(pv[1, 2], Is.EqualTo(0.001_157_4).Within(0.000_000_1));

        pv = Sofa.Plan94(date1, date2, 7); // Uranus
        Assert.That(pv[0, 0], Is.EqualTo(13.210_377).Within(0.000_001));
        Assert.That(pv[0, 1], Is.EqualTo(13.405_065).Within(0.000_001));
        Assert.That(pv[0, 2], Is.EqualTo(5.684_010).Within(0.000_001));
        Assert.That(pv[1, 0], Is.EqualTo(-0.002_934_3).Within(0.000_000_1));
        Assert.That(pv[1, 1], Is.EqualTo(0.002_239_7).Within(0.000_000_1));
        Assert.That(pv[1, 2], Is.EqualTo(0.001_022_5).Within(0.000_000_1));
    }

    #endregion

#region "Astronomy/FundamentalArgs"

// No methods to test.

    #endregion

#region "Astronomy/PrecNutPolar"

// No methods to test. (this may change later)

    #endregion

#region "Astronomy/RotationAndTime"

    [Test]
    public void Era00_Tests() {
        // Era00(): Earth rotation angle (IAU 2000 model).
        // Tests pass: 2023-02-26.
        double date1, date2;
        double era1, era2, era3;
        (date1, date2) = Sofa.Cal2jd(2023, 2, 24);
        era1 = Sofa.Era00(date1, date2);
        Assert.That(era1, Is.EqualTo(2.675_934).Within(0.000_001));

        (date1, date2) = Sofa.Cal2jd(2023, 2, 24);
        date2 += 0.25;
        era2 = Sofa.Era00(date1, date2);
        Assert.That(era2, Is.EqualTo(4.251_031).Within(0.000_001));

        (date1, date2) = Sofa.Cal2jd(2023, 2, 24);
        date2 += 0.75;
        era3 = Sofa.Era00(date1, date2);
        Assert.That(era3, Is.EqualTo(1.118_039).Within(0.000_001));
    }

    [Test]
    public void Gmst06_Tests() {
        // Gmst06(): Greenwich mean sidereal time (IAU 2006 precession).
        // Here using UT1 for TT (error on order of 100 microarcseconds).
        // Tests pass: 2023-02-27.
        double utc1, utc2, ut11, ut12, gmst;
        (utc1, utc2) = Sofa.Cal2jd(2023, 2, 24);
        (ut11, ut12) = Sofa.Utcut1(utc1, utc2, dut1_);
        gmst = Sofa.Gmst06(ut11, ut12, ut11, ut12);
        // Console.WriteLine($"jd:{utc1+utc2}  ut1:{ut11+ut12}  diff:{utc1+utc2-ut11-ut12}");
        // Console.WriteLine($"gmst: {gmst}");
        Assert.That(gmst, Is.EqualTo(2.681_108_9).Within(0.000_000_1));

        ut12 += 0.75;
        gmst = Sofa.Gmst06(ut11, ut12, ut11, ut12);
        // Console.WriteLine();
        // Console.WriteLine($"jd:{utc1+utc2}  ut1:{ut11+ut12}  diff:{utc1+utc2-ut11-ut12}");
        // Console.WriteLine($"gmst: {gmst}");
        Assert.That(gmst, Is.EqualTo(1.123_214_7).Within(0.000_000_1));
    }

    [Test]
    public void Gst06a_Tests() {
        // Gst06a(): Greenwich apparent sidereal time (IAU 2000 and 2006 resolutions).
        // Here using UT1 for TT (error on order of 100 microarcseconds).
        // Tests pass: 2023-02-27.
        double utc1, utc2, ut11, ut12, gast;
        (utc1, utc2) = Sofa.Cal2jd(2023, 2, 24);
        (ut11, ut12) = Sofa.Utcut1(utc1, utc2, dut1_);
        gast = Sofa.Gst06a(ut11, ut12, ut11, ut12);
        Console.WriteLine($"jd:{utc1+utc2}  ut1:{ut11+ut12}  diff:{utc1+utc2-ut11-ut12}");
        Console.WriteLine($"gast: {gast}");
        Assert.That(gast, Is.EqualTo(2.681_068_0).Within(0.000_000_1));

        ut12 += 0.75;
        gast = Sofa.Gst06a(ut11, ut12, ut11, ut12);
        Console.WriteLine();
        Console.WriteLine($"jd:{utc1+utc2}  ut1:{ut11+ut12}  diff:{utc1+utc2-ut11-ut12}");
        Console.WriteLine($"gast: {gast}");
        Assert.That(gast, Is.EqualTo(1.123_173_5).Within(0.000_000_1));
    }

    #endregion

#region "Astronomy/StarCatalogs"

    [Test]
    public void Starpm_Tests() {
        // Starpm(): Update (for new time) star catalog data for proper motion.
        // Tests pass: 2023-02-26.
        var raBefore = HoursAsRadians(ra_yy_gem);
        var decBefore = DegreesAsRadians(dec_yy_gem);
        const double pmr1 = +0.000_001; // radians/year
        const double pmd1 = -0.000_002;
        const double px1 = 0.01;
        const double rv1 = 22.0;
        double utc1, utc2, tai1, tai2, dta, tt1, tt2, ep1a, ep1b, ep2a, ep2b;
        double raAfter, decAfter, pmr2, pmd2, px2, rv2;
        // Construct Before (catalog) epoch:
        (utc1, utc2) = Sofa.Cal2jd(2015, 7, 1);
        (tai1, tai2) = Sofa.Utctai(utc1, utc2);
        dta = Sofa.Dat(2015, 7, 1, 0.0);
        (tt1, tt2) = Sofa.Taiut1(tai1, tai2, dta);
        (ep1a, ep1b) = (tt1, tt2); // omitting dtr (~1 ms, periodic)
        // Construct After (observational) epoch:
        (utc1, utc2) = Sofa.Cal2jd(2023, 2, 24);
        (tai1, tai2) = Sofa.Utctai(utc1, utc2);
        dta = Sofa.Dat(2023, 2, 24, 0.0);
        (tt1, tt2) = Sofa.Taiut1(tai1, tai2, dta);
        (ep2a, ep2b) = (tt1, tt2); // omitting dtr (~1 ms, periodic)
        (raAfter, decAfter, pmr2, pmd2, px2, rv2) = Sofa.Starpm(raBefore, decBefore, pmr1, pmd1, px1, rv1,
            ep1a, ep1b, ep2a, ep2b);
        Assert.That(raBefore, Is.EqualTo(1.983_668_4).Within(0.000_000_1)); // input check
        Assert.That(decBefore, Is.EqualTo(0.556_226_7).Within(0.000_000_1));
        Assert.That(raAfter, Is.EqualTo(1.983_676_1).Within(0.000_000_1));
        Assert.That(decAfter, Is.EqualTo(0.556_211_4).Within(0.000_000_1));

        // Verify exception is thrown whenever Pmsafe() should be used instead:
        var ex = Assert.Throws<ArgumentException>(() => {
            Sofa.Starpm(raBefore, decBefore, pmr1, pmd1, 0.0, 0.0,
                ep1a, ep1b, ep2a, ep2b);
        });
        StringAssert.Contains("Use Sofa.Pmsafe()", ex.Message.ToString());
    }

    #endregion

#region "Astronomy/EclipticCoordinates"

    [Test]
    public void Eceq06_Tests() {
        // Eceq06(): Ecliptic coordinates -> ICRS RA,Dec.
        // Ecliptic coordinate frame is equivalent to astropy's 'barycentricmeanecliptic' with given Time. 
        // Tests pass: 2023-02-26.
        double utc1, utc2, tai1, tai2, dta, tt1, tt2, dr, dd;
        (utc1, utc2) = Sofa.Cal2jd(2023, 2, 24);
        (tai1, tai2) = Sofa.Utctai(utc1, utc2);
        (tt1, tt2) = Sofa.Taiut1(tai1, tai2, dta_);
        var dl = DegreesAsRadians(+111.0);
        var db = DegreesAsRadians(-20.0);
        (dr, dd) = Sofa.Eceq06(tt1, tt2, dl, db);
        Assert.That(dr, Is.EqualTo(DegreesAsRadians(109.391_094_7)).Within(0.000_000_1));
        Assert.That(dd, Is.EqualTo(DegreesAsRadians(2.055_424_7)).Within(0.000_000_1));
    }

    [Test]
    public void Eqec06_Tests() {
        // Eqec06(): ICRS RA,Dec -> Ecliptic coordinates.
        // This ecliptic frame is equivalent to astropy's 'barycentricmeanecliptic' with equinox=Time(). 
        // Tests pass: 2023-02-26.
        double utc1, utc2, tai1, tai2, dta, tt1, tt2, dr, dd, dl, db;
        (utc1, utc2) = Sofa.Cal2jd(2023, 2, 24);
        (tai1, tai2) = Sofa.Utctai(utc1, utc2);
        (tt1, tt2) = Sofa.Taiut1(tai1, tai2, dta_);
        dr = DegreesAsRadians(224.0);
        dd = DegreesAsRadians(55.0);
        (dl, db) = Sofa.Eqec06(tt1, tt2, dr, dd);
        Assert.That(dl, Is.EqualTo(DegreesAsRadians(185.828_904_7)).Within(0.000_000_1));
        Assert.That(db, Is.EqualTo(DegreesAsRadians( 65.511_362_1)).Within(0.000_000_1));
    }
    #endregion

#region "Astronomy/GalacticCoordinates"

[Test]
public void Icrs2g_Tests() {
    // Icrs2g(): ICRS RA,Dec -> Galactic coordinates.
    // Tests pass: 2023-02-26.
    double dr, dd, dl, db;
    dr = DegreesAsRadians(224.0);
    dd = DegreesAsRadians(55.0);
    (dl, db) = Sofa.Icrs2g(dr, dd);
    Assert.That(dl, Is.EqualTo(DegreesAsRadians(92.512_497_6)).Within(0.000_000_1));
    Assert.That(db, Is.EqualTo(DegreesAsRadians(54.138_690_8)).Within(0.000_000_1));
}

[Test]
public void G2icrs_Tests() {
    // G2icrs(): Galactic coordinates -> ICRS RA,Dec.
    // Tests pass: 2023-02-26.
    double dr, dd, dl, db;
    dl = DegreesAsRadians(+111.0);
    db = DegreesAsRadians(-20.0);
    (dr, dd) = Sofa.G2icrs(dl, db);
    Assert.That(dr, Is.EqualTo(DegreesAsRadians(357.830_877_1)).Within(0.000_000_1));
    Assert.That(dd, Is.EqualTo(DegreesAsRadians(41.474_766_4)).Within(0.000_000_1));
}
#endregion
    
#region "Astronomy/GeodeticGeocentric"

[Test]
public void Gc2gd_Tests() {
    // Gc2gd(): Geocentric coordinates (x,y,z) -> Geodetic (long, lat, height),
    //          using specified reference ellipsoid.
    // Tests pass: 2023-02-26.
    int n = 1; // specifies WGS84 reference ellipsoid.
    double[] xyz = new double[3];
    xyz[0] = -0_700_000;
    xyz[1] = +5_580_000;
    xyz[2] = +3_000_000;
    var (elong, phi, height) = Sofa.Gc2gd(n, xyz);
    Assert.That(elong,  Is.EqualTo(DegreesAsRadians(97.150_290)).Within(0.000_010));
    Assert.That(phi,    Is.EqualTo(DegreesAsRadians(28.237_941)).Within(0.000_010));
    Assert.That(height, Is.EqualTo(503.189).Within(0.001));
}

[Test]
public void Gd2gc_Tests() {
    // Gd2gc(): Geodetic (long, lat, height) -> Geocentric coordinates (x,y,z),
    //          using specified reference ellipsoid.
    // Tests pass: 2023-02-26.
    int n = 1; // specifies WGS84reference ellipsoid.
    double elong  = DegreesAsRadians(97);
    double phi    = DegreesAsRadians(28);
    double height = 500; 
    var xyz = new double[3] {0, 0, 0};
    xyz = Sofa.Gd2gc(n, elong, phi, height);
    Assert.That(xyz[0], Is.EqualTo(  -686_875).Within(1));
    Assert.That(xyz[1], Is.EqualTo(+5_594_150).Within(1));
    Assert.That(xyz[2], Is.EqualTo(+2_976_740).Within(1));
}
#endregion

#region "Astronomy/Timescales"

[Test]
public void D2dtf_Tests() {
    // D2dtf(): Format for output a 2-part Julian Date
    // Tests pass: 2023-02-26.
    const string scale = "UTC";
    const int ndp = 6;
    var (d1, d2) = Sofa.Cal2jd(2023, 2, 24);
    d2 += 0.6777;
    // var ihmsf = new int[4] {0, 0, 0, 0};
    var (iy, im, id, ihmsf) = Sofa.D2dtf(scale, ndp, d1, d2);
    Assert.That(iy, Is.EqualTo(2023));
    Assert.That(im, Is.EqualTo(   2));
    Assert.That(id, Is.EqualTo(  24));
    Assert.That(ihmsf[0], Is.EqualTo(16));
    Assert.That(ihmsf[1], Is.EqualTo(15));
    Assert.That(ihmsf[2], Is.EqualTo(53));
    Assert.That(ihmsf[3], Is.EqualTo(280000));
}

[Test]
public void Dat_Tests() {
    // Dat(): Given UTC date, return Delta(AT) = TAI-UTC (cumulative leap seconds).
    // Tests pass: 2023-02-26.
    int iy = 2023, im = 2, id = 24;
    double fd = 0.6777;
    double deltat;
    deltat = Sofa.Dat(iy, im, id, fd);
    Assert.That(deltat, Is.EqualTo(dat_));
    iy = 2003;
    im = 2;
    id = 24;
    deltat = Sofa.Dat(iy, im, id, fd);
    Assert.That(deltat, Is.EqualTo(32));
}

[Test]
public void Dtf2d_Tests() {
    // Dtf2d(): Encode date and time fields into 2−part Julian Date (UTC).
    // Tests pass: 2023-02-26.
    string scale = "UTC";
    int iy = 2023, im = 02, id = 24, ihr = 11, imn = 23;
    double sec = 22.5;
    var (d1, d2) = Sofa.Dtf2d(scale, iy, im, id, ihr, imn, sec);
    Assert.That(d1, Is.EqualTo(2459999.5).Within(0.000_001));
    Assert.That(d2, Is.EqualTo(0.474_566).Within(0.000_001));
}

[Test]
public void Taitt_Tests() {
    // Taitt(): TAI -> TT (adds 32.184 sec).
    // Tests pass: 2023-02-27.
    var (tai1, tai2) = Sofa.Cal2jd(2023, 2, 24);
    var (tt1, tt2) = Sofa.Taitt(tai1, tai2);
    Assert.That(tt1, Is.EqualTo(2_400_000.5).Within(0.000_001));
    Assert.That(tt2, Is.EqualTo(   59_999.000_372_5).Within(0.000_000_1));
}

[Test]
public void Taiut1_Tests() {
    // Taiut1(): TAI -> UT1 (nb: UT1-TAI = "Delta T" = dta)
    // Tests pass: 2023-02-27.
    var (tai1, tai2) = Sofa.Cal2jd(2023, 2, 24);
    var (ut11, ut12) = Sofa.Taiut1(tai1, tai2, dta_);
    Assert.That(ut11, Is.EqualTo(2400000.5).Within(0.000_001));
    Assert.That(ut12, Is.EqualTo(  59998.99957159259).Within(0.000_000_1));
}

[Test]
public void Taiutc_Tests() {
    // Taiut1(): TAI -> UTC (differ by cumulative leap seconds).
    // Tests pass: 2023-02-26.
    var (tai1, tai2) = Sofa.Cal2jd(2023, 2, 24);
    var (utc1, utc2) = Sofa.Taiutc(tai1, tai2);
    Assert.That(tai1, Is.EqualTo(2400000.5).Within(0.000_001));  // input check
    Assert.That(tai2, Is.EqualTo(  59999.0).Within(0.000_001));  // "
    Assert.That(utc1, Is.EqualTo(2400000.5).Within(0.000_001));       // test.
    Assert.That(utc2, Is.EqualTo(  59998.999_572).Within(0.000_001)); // "
}

[Test]
public void Tcbtdb_Tests() {
    // Tcbtdb(): TCB -> TDB ~= Teph (JPL).
    // Tests pass: 2023-02-26 (and matches astropy result).
    var (tcb1, tcb2) = Sofa.Cal2jd(2023, 2, 24);
    var (tdb1, tdb2) = Sofa.Tcbtdb(tcb1, tcb2);
    Assert.That(tcb1, Is.EqualTo(2400000.5).Within(0.000_001));  // input check
    Assert.That(tcb2, Is.EqualTo(  59999.0).Within(0.000_001));  // "
    Assert.That(tdb1, Is.EqualTo(2400000.5).Within(0.000_001));       // test.
    Assert.That(tdb2, Is.EqualTo(  59998.999_739).Within(0.000_001)); // "
}

[Test]
public void Tcgtt_Tests() {
    // Tcgtt(): Geocentric Coordinate Time TCG -> TT.
    // Tests pass: 2023-02-26 (and matches astropy result).
    var (tcg1, tcg2) = Sofa.Cal2jd(2023, 2, 24);
    var (tt1, tt2) = Sofa.Tcgtt(tcg1, tcg2);
    Assert.That(tcg1, Is.EqualTo(2400000.5).Within(0.000_001));  // input check
    Assert.That(tcg2, Is.EqualTo(  59999.0).Within(0.000_001));  // "
    Assert.That(tt1, Is.EqualTo(2400000.5).Within(0.000_001));       // test.
    Assert.That(tt2, Is.EqualTo(  59998.999_988_25).Within(0.000_000_1)); // "
}

[Test]
public void Tdbtcb_Tests() {
    // Tdbtcb(): Barycentric Dynamical Time TDB -> Barycentric Coordinate Time TCB.
    // Tests pass: 2023-02-26 (and matches astropy result).
    var (tdb1, tdb2) = Sofa.Cal2jd(2023, 2, 24);
    var (tcb1, tcb2) = Sofa.Tdbtcb(tdb1, tdb2);
    Assert.That(tdb1, Is.EqualTo(2400000.5).Within(0.000_001));  // input check
    Assert.That(tdb2, Is.EqualTo(  59999.0).Within(0.000_001));  // "
    Assert.That(tcb1, Is.EqualTo(2400000.5).Within(0.000_001));       // test.
    Assert.That(tcb2, Is.EqualTo(  59999.000_261_34).Within(0.000_000_1)); // "
}

[Test]
public void Tttai_Tests() {
    // Tttai(): TT -> TAI (constant difference).
    // Tests pass: 2023-02-26 (and matches astropy result).
    var (tt1, tt2) = Sofa.Cal2jd(2023, 2, 24);
    var (tai1, tai2) = Sofa.Tttai(tt1, tt2);
    Assert.That(tt1,  Is.EqualTo(2400000.5).Within(0.000_001));  // input check
    Assert.That(tt2,  Is.EqualTo(  59999.0).Within(0.000_001));  // "
    Assert.That(tai1, Is.EqualTo(2400000.5).Within(0.000_001));           // test.
    Assert.That(tai2, Is.EqualTo(  59998.999_627_5).Within(0.000_000_1)); // "
}

[Test]
public void Tttcg_Tests() {
    // Tttcg(): TT -> TCG (Geocentric Coordinate Time).
    // Tests pass: 2023-02-27 (and matches astropy result).
    var (tt1, tt2) = Sofa.Cal2jd(2023, 2, 24);
    var (tcg1, tcg2) = Sofa.Tttcg(tt1, tt2);
    Assert.That(tt1,  Is.EqualTo(2400000.5).Within(0.000_001));  // input check
    Assert.That(tt2,  Is.EqualTo(  59999.0).Within(0.000_001));  // "
    Assert.That(tcg1, Is.EqualTo(2400000.5).Within(0.000_001));           // test.
    Assert.That(tcg2, Is.EqualTo(  59999.000_011_7).Within(0.000_000_1)); // "
}

[Test]
public void Ttut1_Tests() {
    // Ttut1(): TT -> UT1, given TT-UT1.
    // Tests pass: 2023-02-27 (and match astropy result).
    var (tt1, tt2) = Sofa.Cal2jd(2023, 2, 24);
    var (ut11, ut12) = Sofa.Ttut1(tt1, tt2, dt_);
    Assert.That(tt1, Is.EqualTo(2400000.5).Within(0.000_001));  // input check
    Assert.That(tt2, Is.EqualTo(59999.0).Within(0.000_001));  // "
    Assert.That(dt_,  Is.EqualTo(69.1984).Within(0.000_1));      // "
    Assert.That(ut11, Is.EqualTo(2400000.5).Within(0.000_001));           // test.
    Assert.That(ut12, Is.EqualTo(59998.999_199_1).Within(0.000_000_1)); // "
}

[Test]
public void Ut1tai_Tests() {
    // Ut1tai(): UT1 -> TAI, given UT1-TAI, available from IERS.
    // Tests pass: 2023-02-27 (and match astropy result).
    var (ut11, ut12) = Sofa.Cal2jd(2023, 2, 24);
    var (tai1, tai2) = Sofa.Ttut1(ut11, ut12, dta_);
    Assert.That(ut11, Is.EqualTo(2400000.5).Within(0.000_001)); // input check
    Assert.That(ut12, Is.EqualTo(59999.0).Within(0.000_001));   // "
    Assert.That(dta_,  Is.EqualTo(-37.0144).Within(0.000_1));     // "
    Assert.That(tai1, Is.EqualTo(2400000.5).Within(0.000_001));         // test.
    Assert.That(tai2, Is.EqualTo(59999.000_428_5).Within(0.000_000_1)); // "
}

[Test]
public void Ut1tt_Tests() {
    // Ut1tt(): UT1 -> TT, given TT-UT1.
    // Tests pass: 2023-02-27 (and match astropy result).
    var (ut11, ut12) = Sofa.Cal2jd(2023, 2, 24);
    var (tt1, tt2) = Sofa.Ut1tt(ut11, ut12, dt_);
    Assert.That(ut11, Is.EqualTo(2400000.5).Within(0.000_001)); // input check
    Assert.That(ut12, Is.EqualTo(59999.0).Within(0.000_001));   // "
    Assert.That(dt_,  Is.EqualTo(69.198_4).Within(0.000_1));     // "
    Assert.That(tt1, Is.EqualTo(2400000.5).Within(0.000_001));         // test.
    Assert.That(tt2, Is.EqualTo(59999.000_800_8).Within(0.000_000_2)); // "
}

[Test]
public void Ut1utc_Tests() {
    // Ut1utc(): UT1 -> UTC (- Earth rotation correction).
    // Tests pass: 2023-02-27 (and match astropy result).
    var (ut11, ut12) = Sofa.Cal2jd(2023, 2, 24);
    var (utc1, utc2) = Sofa.Ut1utc(ut11, ut12, dut1_);
    Assert.That(ut11, Is.EqualTo(2400000.5).Within(0.000_001)); // input check
    Assert.That(ut12, Is.EqualTo(59999.0).Within(0.000_001));   // "
    Assert.That(dut1_,  Is.EqualTo(-0.0144).Within(0.000_1));    // "
    Assert.That(utc1, Is.EqualTo(2400000.5).Within(0.000_001));           // test.
    Assert.That(utc2, Is.EqualTo(59999.000_000_167).Within(0.000_000_1)); // "
}

[Test]
public void Utctai_Tests() {
    // Utctai(): UTC -> TAI (- cumulative leap seconds).
    // Tests pass: 2023-02-27 (and match astropy result).
    var (utc1, utc2) = Sofa.Cal2jd(2023, 2, 24);
    var (tai1, tai2) = Sofa.Utctai(utc1, utc2);
    Assert.That(utc1, Is.EqualTo(2400000.5).Within(0.000_001)); // input check
    Assert.That(utc2, Is.EqualTo(59999.0).Within(0.000_001));   // "
    Assert.That(tai1, Is.EqualTo(2400000.5).Within(0.000_001));         // test.
    Assert.That(tai2, Is.EqualTo(59999.000_428_2).Within(0.000_000_1)); // "
}

[Test]
public void Utcut1_Tests() {
    // Utcut1(): UTC -> UT1 (+ Earth rotation correction).
    // Tests pass: 2023-02-27 (and match astropy result).
    var (utc1, utc2) = Sofa.Cal2jd(2023, 2, 24);
    var (ut11, ut12) = Sofa.Utcut1(utc1, utc2, dut1_);
    Assert.That(utc1, Is.EqualTo(2400000.5).Within(0.000_001)); // input check
    Assert.That(utc2, Is.EqualTo(59999.0).Within(0.000_001));   // "
    Assert.That(dut1_,  Is.EqualTo(-0.0144).Within(0.000_1));    // "
    Assert.That(ut11, Is.EqualTo(2400000.5).Within(0.000_001));           // test.
    Assert.That(ut12, Is.EqualTo(59998.999_999_833).Within(0.000_000_001)); // "
}
#endregion

#region "Astronomy/HorizonEquatorial"

[Test]
public void Ae2hd_Tests() {
    // Ae2hd(): AzAlt -> Hour Angle + (topocentric) Declination.
    // Tests pass: 2023-02-27 (and match TheSkyX result).
    var az = DegreesAsRadians(145.0);
    var el = DegreesAsRadians(67.0); // altitude.
    var phi = DegreesAsRadians(nms_latitude);
    var (ha, dec) = Sofa.Ae2hd(az, el, phi);
    Assert.That(ha, Is.EqualTo(DegreesAsRadians(-0.887_9 * 15.0)).Within(0.000_1));   
    Assert.That(dec, Is.EqualTo(DegreesAsRadians(13.374)).Within(0.000_1));
}

[Test]
public void Hd2ae_Tests() {
    // Hd2ae(): Hour Angle + (topocentric) Declination -> AzAlt.
    // Tests pass: 2023-02-27 (and match TheSkyX result).
    var ha = DegreesAsRadians(15 * (-1.1));  // -1.1 hours
    var dec = DegreesAsRadians(+20);
    var phi = DegreesAsRadians(nms_latitude);
    var (az, el) = Sofa.Hd2ae(ha, dec, phi);
    Assert.That(az, Is.EqualTo(DegreesAsRadians(127.159_51)).Within(0.000_1));   
    Assert.That(el, Is.EqualTo(DegreesAsRadians( 70.434_59)).Within(0.000_1));
}
#endregion

#region "Astronomy/Gnomonic"

[Test]
public void Tpsts_Tests() {
    // Tpsts(): In tangent plane, star xy coord + tangent pt sph coord -> star spherical coordinates.
    // Tests pass: 2023-02-27.
    var xi = DegreesAsRadians(0.22);
    var eta = DegreesAsRadians(0.07); // due North.
    var a0 = DegreesAsRadians(150.0); // RA-like.
    var b0 = DegreesAsRadians(33.0); // Dec-like.
    var (a, b) = Sofa.Tpsts(xi, eta, a0, b0);
    Assert.That(a, Is.EqualTo(DegreesAsRadians(150.262_526)).Within(0.000_001));   
    Assert.That(b, Is.EqualTo(DegreesAsRadians( 33.069_725)).Within(0.000_001));
}

[Test]
public void Tpxes_Tests() {
    // Tpxes(): In tangent plane, star spherical coordinates + tangent pt sph coord -> star xy coord.
    // Tests pass: 2023-02-27.
    var a = DegreesAsRadians(150.0);
    var b = DegreesAsRadians( 61.0);
    var a0 = DegreesAsRadians(149.0); // RA-like.
    var b0 = DegreesAsRadians( 60.0); // Dec-like.
    var (xi, eta) = Sofa.Tpxes(a, b, a0, b0);
    Assert.That(xi,  Is.EqualTo(DegreesAsRadians(0.484_877)).Within(0.000_001));   
    Assert.That(eta, Is.EqualTo(DegreesAsRadians(1.003_803)).Within(0.000_001));
}

#endregion

#region "VectorMatrix/SeparationAndAngle"

[Test]
public void Sepp_Tests() {
    // Sepp(): Angular separation between 2 position (xyz) vectors.
    // Tests pass: 2023-02-27.
    var a = new double[3] {1, 2, 3};
    var b = new double[3] {0, 3, 2};
    var sep = Sofa.Sepp(a, b);
    Assert.That(sep, Is.EqualTo(DegreesAsRadians(27.189_619)).Within(0.000_001));
}

[Test]
public void Seps_Tests() {
    // Seps(): Angular separation between 2 spherical coordinate sets.
    // Tests pass: 2023-02-27.
    var al = DegreesAsRadians(250.5); // longitude
    var ap = DegreesAsRadians(33.0); // latitude
    var bl = DegreesAsRadians(200.0); // longitude
    var bp = DegreesAsRadians(-5.0); // latitude
    var sep = Sofa.Seps(al, ap, bl, bp);
    Assert.That(sep, Is.EqualTo(DegreesAsRadians(61.055_533)).Within(0.000_001));
}

#endregion



// Tests pass: 2023-02-27.

    // Console.WriteLine($"ut1:{ut11} {ut12} -> utc:{utc1} {utc2}");
    // Console.WriteLine($"dut1:{dut1_}");
    // Console.WriteLine($"ut1:{ut11 + ut12} -> utc:{utc1 + utc2}  " +
    // $"diff:{(utc1 + utc2 - ut11 - ut12) * 24 * 3600} seconds.");

// Console.WriteLine($"{pv[0, 0]} {pv[0, 1]} {pv[0, 2]}");
// Console.WriteLine($"{pv[1, 0]} {pv[1, 1]} {pv[1, 2]}");
// Console.WriteLine($"{date1+date2}: era={era2}");
// Console.WriteLine($"diff={(era2-era1)/(2.0 * Math.PI)} rotations.");
// Console.WriteLine($"{utc1+utc2}: " +
// $"Ec.lon.deg={RadiansAsDegrees(dl)} Ec.lat.deg={RadiansAsDegrees(db)} -> " +
//     $"RA.deg={RadiansAsDegrees(dr)} Dec.deg={RadiansAsDegrees(dd)}");

}


