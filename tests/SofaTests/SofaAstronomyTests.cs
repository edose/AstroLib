using NUnit.Framework;
using AstroLib.Core.SOFA;
using AstroLib.Core;
using Microsoft.VisualStudio.TestPlatform.Utilities;

// ReSharper disable IdentifierTypo
// ReSharper disable RedundantExplicitArrayCreation
// ReSharper disable SuggestVarOrType_BuiltInTypes

namespace AstroLibTests.SofaTests; 


[TestFixture]
public class SofaAstronomyTests {
    
    private static double DegreesAsRadians(double deg) { return deg * (Math.PI / 180.0); }

    private static double HoursAsRadians(double hours) { return (hours * 15.0) * (Math.PI / 180.0); }

    private static double RadiansAsDegrees(double radians) { return radians * (180.0 / Math.PI); }

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
    
    // NOTE: Almost all tests here are adapted to AstroLib
    // from http://iausofa.org/2021_0512_C/sofa/t_sofa_c.html.

    #region "Astronomy: Astrometry"
    /*######### Astronomy: Astrometry ####################################################################*/

    [Test]
    public void Ab_Tests() {
        // No involvement of iauAstrom struct.
        var pnat = new double[3] {-0.76321968546737951, -0.60869453983060384, -0.21676408580639883};
        var v = new double[3] {2.1044018893653786e-5, -8.9108923304429319e-5, -3.8633714797716569e-5};
        const double s = 0.99980921395708788;
        const double bm1 = 0.99999999506209258;
        var ppr = Sofa.Ab(pnat, v, s, bm1);
        
        var pprExpected = new double[]
            {-0.7631631094219556269, -0.6087553082505590832, -0.2167926269368471279};
        Assert.That(ppr, Is.EqualTo(pprExpected).Within(1e-12));
    }

    [Test]
    public void Apcg13_Tests() {
        const double date1 = 2456165.5;
        const double date2 = 0.401182685;
        var astrom = Sofa.Apcg13(date1, date2);
        
        Assert.IsInstanceOf(typeof(Sofa.iauASTROM), astrom);
        Assert.That(astrom.pmt, Is.EqualTo(12.65133794027378508).Within(1e-11));
        var ebExpected = new double[] 
            {0.9013108747340644755, -0.4174026640406119957, -0.1809822877867817771};
        Assert.That(astrom.eb, Is.EqualTo(ebExpected).Within(1e-12));
        var ehExpected = new double[] 
            {0.8940025429255499549, -0.4110930268331896318, -0.1782189006019749850};
        Assert.That(astrom.eh, Is.EqualTo(ehExpected).Within(1e-12));
        Assert.That(astrom.em, Is.EqualTo(1.010465295964664178).Within(1e-12));
        var vExpected = new double[] 
            {0.4289638912941341125e-4, 0.8115034032405042132e-4, 0.3517555135536470279e-4};
        Assert.That(astrom.v, Is.EqualTo(vExpected).Within(1e-16));
        Assert.That(astrom.bm1, Is.EqualTo(0.9999999951686013142).Within(1e-12));
        var bpn2Dim = new double[,] {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};
        var bpnExpected = AstroMath.Flatten2DimArray(bpn2Dim);
        Assert.That(astrom.bpn, Is.EqualTo(bpnExpected));
    }

    [Test]
    public void Apci13_Tests() {
        const double date1 = 2456165.5;
        const double date2 = 0.401182685;
        var (astrom, eo) = Sofa.Apci13(date1, date2);
        
        Assert.IsInstanceOf(typeof(Sofa.iauASTROM), astrom);
        Assert.That(astrom.pmt, Is.EqualTo(12.65133794027378508).Within(1e-11));
        var ebExpected = new double[] 
            {0.9013108747340644755, -0.4174026640406119957, -0.1809822877867817771};
        Assert.That(astrom.eb, Is.EqualTo(ebExpected).Within(1e-12));
        var ehExpected = new double[] 
            {0.8940025429255499549, -0.4110930268331896318, -0.1782189006019749850};
        Assert.That(astrom.eh, Is.EqualTo(ehExpected).Within(1e-12));
        Assert.That(astrom.em, Is.EqualTo(1.010465295964664178).Within(1e-12));
        var vExpected = new double[] 
            {0.4289638912941341125e-4, 0.8115034032405042132e-4, 0.3517555135536470279e-4};
        Assert.That(astrom.v, Is.EqualTo(vExpected).Within(1e-16));
        Assert.That(astrom.bm1, Is.EqualTo(0.9999999951686013142).Within(1e-12));
        // Transposed here, because SOFA C++ test code gives 2-D array in non-C++ (Fortran) order.
        var bpn2Dim = AstroMath.TransposeSquareMatrix(new double[,] 
            {{0.9999992060376761710, 0.4124244860106037157e-7, 0.1260128571051709670e-2}, 
             {-0.1282291987222130690e-7, 0.9999999997456835325, -0.2255288829420524935e-4}, 
             {-0.1260128571661374559e-2, 0.2255285422953395494e-4, 0.9999992057833604343}});
        var bpnExpected = AstroMath.Flatten2DimArray(bpn2Dim);
        Assert.That(astrom.bpn, Is.EqualTo(bpnExpected).Within(1e-12));
        Assert.That(eo, Is.EqualTo(-0.2900618712657375647e-2).Within(1e-12));
    }

    [Test]
    public void Apco13_Tests() {
        double utc1 = 2456384.5;
        double utc2 = 0.969254051;
        double dut1 = 0.1550675;
        double elong = -0.527800806;
        double phi = -1.2345856;
        double hm = 2738.0;
        double xp = 2.47230737e-7;
        double yp = 1.82640464e-6;
        double phpa = 731.0;
        double tc = 12.8;
        double rh = 0.59;
        double wl = 0.55;
        var (status, astrom, eo) = 
            Sofa.Apco13(utc1, utc2, dut1, elong, phi, hm, xp, yp, phpa, tc, rh, wl);
        
        Assert.IsInstanceOf(typeof(Sofa.iauASTROM), astrom);
        Assert.That(status, Is.EqualTo(0));
        Assert.That(astrom.pmt, Is.EqualTo(13.25248468622475727).Within(1e-11));
        var ebExpected = new double[] 
            {-0.9741827107320875162, -0.2115130190489716682, -0.09179840189496755339};
        Assert.That(astrom.eb, Is.EqualTo(ebExpected).Within(1e-12));
        var ehExpected = new double[] 
            {-0.9736425572586935247, -0.2092452121603336166, -0.09075578153885665295};
        Assert.That(astrom.eh, Is.EqualTo(ehExpected).Within(1e-12));
        Assert.That(astrom.em, Is.EqualTo(0.9998233240913898141).Within(1e-12));
        var vExpected = new double[] 
            {0.2078704994520489246e-4, -0.8955360133238868938e-4, -0.3863338993055887398e-4};
        Assert.That(astrom.v, Is.EqualTo(vExpected).Within(1e-16));
        Assert.That(astrom.bm1, Is.EqualTo(0.9999999950277561004).Within(1e-12));
        // Transposed here, because SOFA C++ test code gives 2-D array in non-C++ (Fortran) order.
        var bpn2Dim = AstroMath.TransposeSquareMatrix(new double[,] 
        {{0.9999991390295147999, 0.4978650075315529277e-7, 0.001312227200850293372}, 
            {-0.1136336652812486604e-7, 0.9999999995713154865, -0.2928086230975367296e-4}, 
            {-0.001312227201745553566, 0.2928082218847679162e-4, 0.9999991386008312212}});
        var bpnExpected = AstroMath.Flatten2DimArray(bpn2Dim);
        Assert.That(astrom.bpn,    Is.EqualTo(bpnExpected).Within(1e-12));
        Assert.That(astrom.along,  Is.EqualTo(-0.5278008060295995733).Within(1e-12));
        Assert.That(astrom.xpl,    Is.EqualTo(0.1133427418130752958e-5).Within(1e-17));
        Assert.That(astrom.ypl,    Is.EqualTo(0.1453347595780646207e-5).Within(1e-17));
        Assert.That(astrom.sphi,   Is.EqualTo(-0.9440115679003211329).Within(1e-12));
        Assert.That(astrom.cphi,   Is.EqualTo(0.3299123514971474711).Within(1e-12));
        Assert.That(astrom.diurab, Is.EqualTo(0.0));
        Assert.That(astrom.eral,   Is.EqualTo(2.617608909189664000).Within(1e-12));
        Assert.That(astrom.refa,   Is.EqualTo(0.2014187785940396921e-3).Within(1e-15));
        Assert.That(astrom.refb,   Is.EqualTo(-0.2361408314943696227e-6).Within(1e-12));
        Assert.That(eo,            Is.EqualTo(-0.003020548354802412839).Within(1e-14));
    }

    [Test]
    public void Apcs13_Tests() {
        double date1 = 2456165.5, date2 = 0.401182685;
        var pv = new double[2, 3] {{-6241497.16, 401346.896, -1251136.04},
                                   {-29.264597, -455.021831, 0.0266151194}};
        var astrom = Sofa.Apcs13(date1, date2, pv);
        
        Assert.IsInstanceOf(typeof(Sofa.iauASTROM), astrom);
        Assert.That(astrom.pmt, Is.EqualTo(12.65133794027378508).Within(1e-11));
        var ebExpected = new double[] 
            {0.9012691529025250644, -0.4173999812023194317, -0.1809906511146429670};
        Assert.That(astrom.eb, Is.EqualTo(ebExpected).Within(1e-12));
        var ehExpected = new double[] 
            {0.8939939101760130792, -0.4111053891734021478, -0.1782336880636997374};
        Assert.That(astrom.eh, Is.EqualTo(ehExpected).Within(1e-12));
        Assert.That(astrom.em, Is.EqualTo(1.010428384373491095).Within(1e-12));
        var vExpected = new double[] 
            {0.4279877294121697570e-4, 0.7963255087052120678e-4, 0.3517564013384691531e-4};
        Assert.That(astrom.v, Is.EqualTo(vExpected).Within(1e-16));
        Assert.That(astrom.bm1, Is.EqualTo(0.9999999952947980978).Within(1e-12));
        var bpn2Dim = AstroMath.TransposeSquareMatrix(new double[,] 
            {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}});
        var bpnExpected = AstroMath.Flatten2DimArray(bpn2Dim);
        Assert.That(astrom.bpn, Is.EqualTo(bpnExpected));
    }

    [Test]
    public void Aper13_Tests() {
        Sofa.iauASTROM astrom = default;  // simulating a pre-existing iauASTROM instance.
        astrom.along = 1.234;
        double ut11 = 2456165.5, ut12 = 0.401182685;
        Assert.That(astrom.eral, Is.EqualTo(0));  // because test struct is newly initialized.
        Sofa.Aper13(ut11, ut12, ref astrom);
        
        Assert.IsInstanceOf(typeof(Sofa.iauASTROM), astrom);
        Assert.That(astrom.eral, Is.EqualTo(3.316236661789694933).Within(1e-12));  // verify the update.
    }

    [Test]
    public void Apio13_Tests() {
        double utc1 = 2456384.5;
        double utc2 = 0.969254051;
        double dut1 = 0.1550675;
        double elong = -0.527800806;
        double phi = -1.2345856;
        double hm = 2738.0;
        double xp = 2.47230737e-7;
        double yp = 1.82640464e-6;
        double phpa = 731.0;
        double tc = 12.8;
        double rh = 0.59;
        double wl = 0.55;
        var (status, astrom) = Sofa.Apio13(utc1, utc2, dut1, elong, phi, hm, xp, yp,
            phpa, tc, rh, wl);
        
        Assert.That(status, Is.EqualTo(0));
        Assert.IsInstanceOf(typeof(Sofa.iauASTROM), astrom);
        Assert.That(astrom.along,  Is.EqualTo(-0.5278008060295995733).Within(1e-12));
        Assert.That(astrom.xpl,    Is.EqualTo(0.1133427418130752958e-5).Within(1e-17));
        Assert.That(astrom.ypl,    Is.EqualTo(0.1453347595780646207e-5).Within(1e-17));
        Assert.That(astrom.sphi,   Is.EqualTo(-0.9440115679003211329).Within(1e-12));
        Assert.That(astrom.cphi,   Is.EqualTo(0.3299123514971474711).Within(1e-12));
        Assert.That(astrom.diurab, Is.EqualTo(0.5135843661699913529e-6).Within(1e-12));
        Assert.That(astrom.eral,   Is.EqualTo(2.617608909189664000).Within(1e-12));
        Assert.That(astrom.refa,   Is.EqualTo(0.2014187785940396921e-3).Within(1e-15));
        Assert.That(astrom.refb,   Is.EqualTo(-0.2361408314943696227e-6).Within(1e-18));
    }

    [Test]
    public void Atcc13_Tests() {
        double rc = 2.71;
        double dc = 0.174;
        double pr = 1e-5;
        double pd = 5e-6;
        double px = 0.1;
        double rv = 55.0;
        double date1 = 2456165.5;
        double date2 = 0.401182685;
        var (ra, da) = Sofa.Atcc13(rc, dc, pr, pd, px, rv, date1, date2);
        
        Assert.That(ra,Is.EqualTo(2.710126504531372384).Within(1e-12));
        Assert.That(da,Is.EqualTo(0.1740632537628350152).Within(1e-12));
    }
    
    [Test]
    public void Atccq_Tests() {
        double date1 = 2456165.5, date2 = 0.401182685;
        var (astrom, eo) = Sofa.Apci13(date1, date2);
        double rc = 2.71;
        double dc = 0.174;
        double pr = 1e-5;
        double pd = 5e-6;
        double px = 0.1;
        double rv = 55.0;
        var (ra, da) = Sofa.Atccq(rc, dc, pr, pd, px, rv, astrom);
        
        Assert.That(ra,Is.EqualTo(2.710126504531372384).Within(1e-12));
        Assert.That(da,Is.EqualTo(0.1740632537628350152).Within(1e-12));
    }

    [Test]
    public void Atci13_Tests() {
        double rc = 2.71;
        double dc = 0.174;
        double pr = 1e-5;
        double pd = 5e-6;
        double px = 0.1;
        double rv = 55.0;
        double date1 = 2456165.5;
        double date2 = 0.401182685;
        var (ri, di, eo) = Sofa.Atci13(rc, dc, pr, pd, px, rv, date1, date2);
        
        Assert.That(ri,Is.EqualTo(2.710121572968696744).Within(1e-12));
        Assert.That(di,Is.EqualTo(0.1729371367219539137).Within(1e-12));
        Assert.That(eo,Is.EqualTo(-0.002900618712657375647).Within(1e-14));
    }

    [Test]
    public void Atciq_Tests() {
        double date1 = 2456165.5, date2 = 0.401182685;
        var (astrom, eo) = Sofa.Apci13(date1, date2);
        double rc = 2.71;
        double dc = 0.174;
        double pr = 1e-5;
        double pd = 5e-6;
        double px = 0.1;
        double rv = 55.0;
        var (ri, di) = Sofa.Atciq(rc, dc, pr, pd, px, rv, astrom);
        
        Assert.That(ri,Is.EqualTo(2.710121572968696744).Within(1e-12));
        Assert.That(di,Is.EqualTo(0.1729371367219539137).Within(1e-12));
    }

    public void Atciqz_Tests() {
        double date1 = 2456165.5, date2 = 0.401182685;
        var (astrom, eo) = Sofa.Apci13(date1, date2);
        double rc = 2.71;
        double dc = 0.174;
        var (ri, di) = Sofa.Atciqz(rc, dc, astrom);
        
        Assert.That(ri,Is.EqualTo(2.709994899247256984).Within(1e-12));
        Assert.That(di,Is.EqualTo(0.1728740720984931891).Within(1e-12));
    }

    [Test]
    public void Atco13_Tests() {
        double rc = 2.71;
        double dc = 0.174;
        double pr = 1e-5;
        double pd = 5e-6;
        double px = 0.1;
        double rv = 55.0;
        double utc1 = 2456384.5;
        double utc2 = 0.969254051;
        double dut1 = 0.1550675;
        double elong = -0.527800806;
        double phi = -1.2345856;
        double hm = 2738.0;
        double xp = 2.47230737e-7;
        double yp = 1.82640464e-6;
        double phpa = 731.0;
        double tc = 12.8;
        double rh = 0.59;
        double wl = 0.55;
        var (status, aob, zob, hob, dob, rob, eo) =
            Sofa.Atco13(rc, dc, pr, pd, px, rv, utc1, utc2, dut1, elong, phi, hm, xp, yp, phpa, tc, rh, wl);
        
        Assert.That(status,Is.EqualTo(0));
        Assert.That(aob,Is.EqualTo(0.9251774485485515207e-1).Within(1e-12));
        Assert.That(zob,Is.EqualTo(1.407661405256499357).Within(1e-12));
        Assert.That(hob,Is.EqualTo(-0.9265154431529724692e-1).Within(1e-12));
        Assert.That(dob,Is.EqualTo(0.1716626560072526200).Within(1e-12));
        Assert.That(rob,Is.EqualTo(2.710260453504961012).Within(1e-12));
        Assert.That(eo,Is.EqualTo(-0.003020548354802412839).Within(1e-14));
    }

    [Test]
    public void Atic13_Tests() {
        double ri = 2.710121572969038991;
        double di = 0.1729371367218230438;
        double date1 = 2456165.5;
        double date2 = 0.401182685;
        var (rc, dc, eo) = Sofa.Atic13(ri, di, date1, date2);

        Assert.That(rc,Is.EqualTo(2.710126504531716819).Within(1e-12));
        Assert.That(dc,Is.EqualTo(0.1740632537627034482).Within(1e-12));
        Assert.That(eo,Is.EqualTo(-0.002900618712657375647).Within(1e-14));
    }

    [Test]
    public void Aticq_Tests() {
        double date1 = 2456165.5, date2 = 0.401182685;
        var (astrom, eo) = Sofa.Apci13(date1, date2);
        double ri = 2.710121572969038991;
        double di = 0.1729371367218230438;
        var (rc, dc) = Sofa.Aticq(ri, di, astrom);
        
        Assert.That(rc,Is.EqualTo(2.710126504531716819).Within(1e-12));
        Assert.That(dc,Is.EqualTo(0.1740632537627034482).Within(1e-12));
    }

    [Test]
    public void Atio13_Tests() {
        double ri = 2.710121572969038991;
        double di = 0.1729371367218230438;
        double utc1 = 2456384.5;
        double utc2 = 0.969254051;
        double dut1 = 0.1550675;
        double elong = -0.527800806;
        double phi = -1.2345856;
        double hm = 2738.0;
        double xp = 2.47230737e-7;
        double yp = 1.82640464e-6;
        double phpa = 731.0;
        double tc = 12.8;
        double rh = 0.59;
        double wl = 0.55;
        var (status, aob, zob, hob, dob, rob) =
            Sofa.Atio13(ri, di, utc1, utc2, dut1, elong, phi, hm, xp, yp, phpa, tc, rh, wl);
        
        Assert.That(status,Is.EqualTo(0));
        Assert.That(aob,Is.EqualTo(0.9233952224895122499e-1).Within(1e-12));
        Assert.That(zob,Is.EqualTo(1.407758704513549991).Within(1e-12));
        Assert.That(hob,Is.EqualTo(-0.9247619879881698140e-1).Within(1e-12));
        Assert.That(dob,Is.EqualTo(0.1717653435756234676).Within(1e-12));
        Assert.That(rob,Is.EqualTo(2.710085107988480746).Within(1e-12));
    }

    [Test]
    public void Atioq_Tests() {
        double utc1 = 2456384.5;
        double utc2 = 0.969254051;
        double dut1 = 0.1550675;
        double elong = -0.527800806;
        double phi = -1.2345856;
        double hm = 2738.0;
        double xp = 2.47230737e-7;
        double yp = 1.82640464e-6;
        double phpa = 731.0;
        double tc = 12.8;
        double rh = 0.59;
        double wl = 0.55;
        var (status, astrom) = Sofa.Apio13(utc1, utc2, dut1, elong, phi, hm, xp, yp, phpa, tc, rh, wl);
        double ri = 2.710121572969038991;
        double di = 0.1729371367218230438;
        var (aob, zob, hob, dob, rob) = Sofa.Atioq(ri, di, astrom);
        
        Assert.That(aob,Is.EqualTo(0.9233952224895122499e-1).Within(1e-12));
        Assert.That(zob,Is.EqualTo(1.407758704513549991).Within(1e-12));
        Assert.That(hob,Is.EqualTo(-0.9247619879881698140e-1).Within(1e-12));
        Assert.That(dob,Is.EqualTo(0.1717653435756234676).Within(1e-12));
        Assert.That(rob,Is.EqualTo(2.710085107988480746).Within(1e-12));
    }

    [Test]
    public void Atoc13_Tests() {
        // Test with coordinate type = "R" (Right ascension, declination):
        double ob1 = 2.710085107986886201;
        double ob2 = 0.1717653435758265198;
        double utc1 = 2456384.5;
        double utc2 = 0.969254051;
        double dut1 = 0.1550675;
        double elong = -0.527800806;
        double phi = -1.2345856;
        double hm = 2738.0;
        double xp = 2.47230737e-7;
        double yp = 1.82640464e-6;
        double phpa = 731.0;
        double tc = 12.8;
        double rh = 0.59;
        double wl = 0.55;
        var (status, rc, dc) =
            Sofa.Atoc13("R", ob1, ob2, utc1, utc2, dut1, elong, phi, hm, xp, yp, phpa, tc, rh, wl);
        Assert.That(status,Is.EqualTo(0));
        Assert.That(rc,Is.EqualTo(2.709956744659136129).Within(1e-12));
        Assert.That(dc,Is.EqualTo(0.1741696500898471362).Within(1e-12));
        
        // Test with coordinate type = "H" (Hour angle, declination):
        ob1 = -0.09247619879782006106;
        ob2 = 0.1717653435758265198;
        (status, rc, dc) =
            Sofa.Atoc13("H", ob1, ob2, utc1, utc2, dut1, elong, phi, hm, xp, yp, phpa, tc, rh, wl);
        Assert.That(status,Is.EqualTo(0));
        Assert.That(rc,Is.EqualTo(2.709956744659734086).Within(1e-12));
        Assert.That(dc,Is.EqualTo(0.1741696500898471362).Within(1e-12));
        
        // Test with coordinate type = "A" (Azimuth, zenith distance):
        ob1 = 0.09233952224794989993;
        ob2 = 1.407758704513722461;
        (status, rc, dc) =
            Sofa.Atoc13("A", ob1, ob2, utc1, utc2, dut1, elong, phi, hm, xp, yp, phpa, tc, rh, wl);
        Assert.That(status,Is.EqualTo(0));
        Assert.That(rc,Is.EqualTo(2.709956744659734086).Within(1e-12));
        Assert.That(dc,Is.EqualTo(0.1741696500898471366).Within(1e-12));
    }

    [Test]
    public void Atoi13_Tests() {
        double utc1, utc2, dut1, elong, phi, hm, xp, yp, phpa, tc, rh, wl, ob1, ob2;
        
        // Test with coordinate type = "R" (Right ascension, declination):
        utc1 = 2456384.5;
        utc2 = 0.969254051;
        dut1 = 0.1550675;
        elong = -0.527800806;
        phi = -1.2345856;
        hm = 2738.0;
        xp = 2.47230737e-7;
        yp = 1.82640464e-6;
        phpa = 731.0;
        tc = 12.8;
        rh = 0.59;
        wl = 0.55;
        ob1 = 2.710085107986886201;
        ob2 = 0.1717653435758265198;
        var (status, ri, di) =
            Sofa.Atoi13("R", ob1, ob2, utc1, utc2, dut1, elong, phi, hm, xp, yp, phpa, tc, rh, wl);
        Assert.That(status,Is.EqualTo(0));
        Assert.That(ri,Is.EqualTo(2.710121574447540810).Within(1e-12));
        Assert.That(di,Is.EqualTo(0.1729371839116608778).Within(1e-12));
        
        // Test with coordinate type = "H" (Hour angle, declination):
        ob1 = -0.09247619879782006106;
        ob2 = 0.1717653435758265198;
        (status, ri, di) =
            Sofa.Atoi13("H", ob1, ob2, utc1, utc2, dut1, elong, phi, hm, xp, yp, phpa, tc, rh, wl);
        Assert.That(status,Is.EqualTo(0));
        Assert.That(ri,Is.EqualTo(2.710121574448138676).Within(1e-12));
        Assert.That(di,Is.EqualTo(0.1729371839116608778).Within(1e-12));
        
        // Test with coordinate type = "A" (Azimuth, zenith distance):
        ob1 = 0.09233952224794989993;
        ob2 = 1.407758704513722461;
        (status, ri, di) =
            Sofa.Atoi13("A", ob1, ob2, utc1, utc2, dut1, elong, phi, hm, xp, yp, phpa, tc, rh, wl);
        Assert.That(status,Is.EqualTo(0));
        Assert.That(ri,Is.EqualTo(2.710121574448138676).Within(1e-12));
        Assert.That(di,Is.EqualTo(0.1729371839116608781).Within(1e-12));
    }

    [Test]
    public void Atoiq_Tests() {
        double utc1, utc2, dut1, elong, phi, hm, xp, yp, phpa, tc, rh, wl;
        utc1 = 2456384.5;
        utc2 = 0.969254051;
        dut1 = 0.1550675;
        elong = -0.527800806;
        phi = -1.2345856;
        hm = 2738.0;
        xp = 2.47230737e-7;
        yp = 1.82640464e-6;
        phpa = 731.0;
        tc = 12.8;
        rh = 0.59;
        wl = 0.55;
        var (status, astrom) = Sofa.Apio13(utc1, utc2, dut1, elong, phi, hm, xp, yp, phpa, tc, rh, wl);
        
        // Test with coordinate type = "R" (Right ascension, declination):
        double ob1 = 2.710085107986886201;
        double ob2 = 0.1717653435758265198;
        var (ri, di) = Sofa.Atoiq("R", ob1, ob2, astrom);
        Assert.That(ri,Is.EqualTo(2.710121574447540810).Within(1e-12));
        Assert.That(di,Is.EqualTo(0.17293718391166087785).Within(1e-12));
        
        // Test with coordinate type = "H" (Hour angle, declination):
        ob1 = -0.09247619879782006106;
        ob2 = 0.1717653435758265198;
        (ri, di) = Sofa.Atoiq("H", ob1, ob2, astrom);
        Assert.That(ri,Is.EqualTo(2.710121574448138676).Within(1e-12));
        Assert.That(di,Is.EqualTo(0.1729371839116608778).Within(1e-12));
        
        // Test with coordinate type = "A" (Azimuth, zenith distance):
        ob1 = 0.09233952224794989993;
        ob2 = 1.407758704513722461;
        (ri, di) = Sofa.Atoiq("A", ob1, ob2, astrom);
        Assert.That(ri,Is.EqualTo(2.710121574448138676).Within(1e-12));
        Assert.That(di,Is.EqualTo(0.1729371839116608781).Within(1e-12));
    }

    [Test]
    public void Pmsafe_Tests() {
        double ra1, dec1, pmr1, pmd1, px1, rv1, ep1a, ep1b, ep2a, ep2b;
        ra1 = 1.234;
        dec1 = 0.789;
        pmr1 = 1e-5;
        pmd1 = -2e-5;
        px1 = 1e-2;
        rv1 = 10.0;
        ep1a = 2400000.5;
        ep1b = 48348.5625;
        ep2a = 2400000.5;
        ep2b = 51544.5;
        var (status, ra2, dec2, pmr2, pmd2, px2, rv2) =
            Sofa.Pmsafe(ra1, dec1, pmr1, pmd1, px1, rv1, ep1a, ep1b, ep2a, ep2b);
        Assert.That(status,Is.EqualTo(0));
        Assert.That(ra2,  Is.EqualTo(1.234087484501017061).Within(1e-12));
        Assert.That(dec2, Is.EqualTo(0.7888249982450468567).Within(1e-12));
        Assert.That(pmr2, Is.EqualTo(0.9996457663586073988e-5).Within(1e-12));
        Assert.That(pmd2, Is.EqualTo(-0.2000040085106754565e-4).Within(1e-16));
        Assert.That(px2,  Is.EqualTo(0.9999997295356830666e-2).Within(1e-12));
        Assert.That(rv2,  Is.EqualTo(10.38468380293920069).Within(1e-10));
    }

    [Test]
    public void Pvtob_Tests() {
        double elong, phi, hm, xp, yp, sp, theta;
        elong = 2.0;
        phi = 0.5;
        hm = 3000.0;
        xp = 1e-6;
        yp = -0.5e-6;
        sp = 1e-8;
        theta = 5.0;
        var pv = Sofa.Pvtob(elong, phi, hm, xp, yp, sp, theta);
        Assert.That(pv.GetLength(0), Is.EqualTo(2));
        Assert.That(pv.GetLength(1), Is.EqualTo(3));
        Assert.That(pv[0,0], Is.EqualTo(4225081.367071159207).Within(1e-5));
        Assert.That(pv[0,1], Is.EqualTo(3681943.215856198144).Within(1e-5));
        Assert.That(pv[0,2], Is.EqualTo(3041149.399241260785).Within(1e-5));
        Assert.That(pv[1,0], Is.EqualTo(-268.4915389365998787).Within(1e-9));
        Assert.That(pv[1,1], Is.EqualTo(308.0977983288903123).Within(1e-9));
        Assert.That(pv[1,2], Is.EqualTo(0.0));
    }

    [Test]
    public void Refco_Tests() {
        double phpa, tc, rh, wl;
        phpa = 800.0;
        tc = 10.0;
        rh = 0.9;
        wl = 0.4;
        var (refa, refb) = Sofa.Refco(phpa, tc, rh, wl);
        
        Assert.That(refa,  Is.EqualTo(0.2264949956241415009e-3).Within(1e-15));
        Assert.That(refb,  Is.EqualTo(-0.2598658261729343970e-6).Within(1e-18));
    }

    public void Starpm_Tests() {
        double ra1, dec1, pmr1, pmd1, px1, rv1, ep1a, ep1b, ep2a, ep2b;
        ra1  =   0.01686756;
        dec1 = -1.093989828;
        pmr1 = -1.78323516e-5;
        pmd1 =  2.336024047e-6;
        px1  =   0.74723;
        rv1  = -21.6;
        ep1a = 2400000.5;
        ep1b = 50083.0;
        ep2a = 2400000.5;
        ep2b = 53736.0;
        var (status, ra2, dec2, pmr2, pmd2, px2, rv2) =
            Sofa.Starpm(ra1, dec1, pmr1, pmd1, px1, rv1, ep1a, ep1b, ep2a, ep2b);
        
        Assert.That(ra2,  Is.EqualTo(0.01668919069414256149).Within(1e-13));
        Assert.That(dec2, Is.EqualTo(-1.093966454217127897).Within(1e-13));
        Assert.That(pmr2, Is.EqualTo(-0.1783662682153176524e-4).Within(1e-17));
        Assert.That(pmd2, Is.EqualTo(0.2338092915983989595e-5).Within(1e-17));
        Assert.That(px2,  Is.EqualTo(0.7473533835317719243).Within(1e-13));
        Assert.That(rv2,  Is.EqualTo(-21.59905170476417175).Within(1e-11));
    }

    #endregion
    
    #region "Astronomy: Calendars"
    /*######### Astronomy: Calendars #####################################################################*/

    [Test]
    public void Cal2jd_Tests() {
        var (status, djm0, djm) = Sofa.Cal2jd(2003, 06, 01);
        Assert.That(status, Is.EqualTo(0));
        Assert.That(djm0,   Is.EqualTo(2400000.5));
        Assert.That(djm,    Is.EqualTo(52791.0));
    }

    [Test]
    public void Jd2cal_Tests() {
        double dj1 = 2400000.5, dj2 = 50123.9999;
        var (status, iy, im, id, fd) = Sofa.Jd2cal(dj1, dj2);
        Assert.That(status, Is.EqualTo(0));
        Assert.That(iy,   Is.EqualTo(1996));
        Assert.That(im,   Is.EqualTo(2));
        Assert.That(id,   Is.EqualTo(10));
        Assert.That(fd,   Is.EqualTo(0.9999).Within(1e-7));
    }

    [Test]
    public void Jdcalf_Tests() {
        double dj1 = 2400000.5, dj2 = 50123.9999;
        var (status, iymdf) = Sofa.Jdcalf(4, dj1, dj2);
        
        Assert.That(status, Is.EqualTo(0));
        var iymdfExpected = new int[] {1996, 2, 10, 9999};
        Assert.That(iymdf,   Is.EqualTo(iymdfExpected));
    }
    
    #endregion
    
    #region "Astronomy: Time Scales"
    /*######### Astronomy: Time Scales ###################################################################*/

    [Test]
    public void D2dtf_Tests() {
        var (status, iy, im, id, ihmsf) = 
            Sofa.D2dtf("UTC", 5, 2400000.5, 49533.99999);
        Assert.That(status, Is.EqualTo(0));
        Assert.That(iy,   Is.EqualTo(1994));
        Assert.That(im,   Is.EqualTo(6));
        Assert.That(id,   Is.EqualTo(30));
        var ihmsfExpected = new int[] {23, 59, 60, 13599};
        Assert.That(ihmsf,   Is.EqualTo(ihmsfExpected));
    }

    [Test]
    public void Dat_Tests() {
        var (status, deltat) = Sofa.Dat(2003, 6, 1, 0.0);
        Assert.That(status, Is.EqualTo(0));
        Assert.That(deltat, Is.EqualTo(32.0));
        
        (status, deltat) = Sofa.Dat(2008, 1, 17, 0.0);
        Assert.That(status, Is.EqualTo(0));
        Assert.That(deltat, Is.EqualTo(33.0));

        (status, deltat) = Sofa.Dat(2017, 9, 1, 0.0);
        Assert.That(status, Is.EqualTo(0));
        Assert.That(deltat, Is.EqualTo(37.0));
    }

    [Test]
    public void Dtdb_Tests() {
        var diff = Sofa.Dtdb(2448939.5, 0.123, 0.76543, 5.0123, 5525.242, 3190.0);
        Assert.That(diff, Is.EqualTo(-0.1280368005936998991e-2).Within(1e-15));
    }

    [Test]
    public void Dtf2d_Tests() {
        var (status, d1, d2) = 
            Sofa.Dtf2d("UTC", 1994, 6, 30, 23, 59, 60.13599);
        Assert.That(status, Is.EqualTo(0));
        Assert.That(d1+d2, Is.EqualTo(2449534.49999).Within(1e-6));
    }

    [Test]
    public void Taitt_Tests() {
        var (status, tt1, tt2) = Sofa.Taitt(2453750.5, 0.892482639);
        Assert.That(tt1, Is.EqualTo(2453750.5).Within(1e-6));
        Assert.That(tt2, Is.EqualTo(0.892855139).Within(1e-12));
    }

    [Test]
    public void Taiut1_Tests() {
        var (status, ut11, ut12) = Sofa.Taiut1(2453750.5, 0.892482639, -32.6659);
        Assert.That(ut11, Is.EqualTo(2453750.5).Within(1e-6));
        Assert.That(ut12, Is.EqualTo(0.8921045614537037037).Within(1e-12));
    }

    [Test]
    public void Taiutc_Tests() {
        var (status, utc1, utc2) = Sofa.Taiutc(2453750.5, 0.892482639);
        Assert.That(utc1, Is.EqualTo(2453750.5).Within(1e-6));
        Assert.That(utc2, Is.EqualTo(0.8921006945555555556).Within(1e-12));
    }

    [Test]
    public void Tcbtdb_Tests() {
        var (status, tdb1, tdb2) = Sofa.Tcbtdb(2453750.5, 0.893019599);
        Assert.That(tdb1, Is.EqualTo(2453750.5).Within(1e-6));
        Assert.That(tdb2, Is.EqualTo(0.8928551362746343397).Within(1e-12));
    }
        
    // TODO: Resume here...
    
    [Test]
    public void Tcgtt_Tests() { }
        
    [Test]
    public void Tdbtcb_Tests() { }
        
    [Test]
    public void Tdbtt_Tests() { }
        
    [Test]
    public void Tttai_Tests() { }
        
    [Test]
    public void Tttcg_Tests() { }
            
    [Test]
    public void Tttdb_Tests() { }
            
    [Test]
    public void Ttut1_Tests() { }
            
    [Test]
    public void Ut1tai_Tests() { }
            
    [Test]
    public void Ut1tt_Tests() { }
            
    [Test]
    public void Ut1utc_Tests() { }
            
    [Test]
    public void Utctai_Tests() { }
            
    [Test]
    public void Utcut1_Tests() { }
    
    
    #endregion
    

}