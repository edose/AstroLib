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
        Assert.That(status, Is.EqualTo(0));
        Assert.That(tt1, Is.EqualTo(2453750.5).Within(1e-6));
        Assert.That(tt2, Is.EqualTo(0.892855139).Within(1e-12));
    }

    [Test]
    public void Taiut1_Tests() {
        var (status, ut11, ut12) = Sofa.Taiut1(2453750.5, 0.892482639, -32.6659);
        Assert.That(status, Is.EqualTo(0));
        Assert.That(ut11, Is.EqualTo(2453750.5).Within(1e-6));
        Assert.That(ut12, Is.EqualTo(0.8921045614537037037).Within(1e-12));
    }

    [Test]
    public void Taiutc_Tests() {
        var (status, utc1, utc2) = Sofa.Taiutc(2453750.5, 0.892482639);
        Assert.That(status, Is.EqualTo(0));
        Assert.That(utc1, Is.EqualTo(2453750.5).Within(1e-6));
        Assert.That(utc2, Is.EqualTo(0.8921006945555555556).Within(1e-12));
    }

    [Test]
    public void Tcbtdb_Tests() {
        var (status, tdb1, tdb2) = Sofa.Tcbtdb(2453750.5, 0.893019599);
        Assert.That(status, Is.EqualTo(0));
        Assert.That(tdb1, Is.EqualTo(2453750.5).Within(1e-6));
        Assert.That(tdb2, Is.EqualTo(0.8928551362746343397).Within(1e-12));
    }

    [Test]
    public void Tcgtt_Tests() {
        var (status, tt1, tt2) = Sofa.Tcgtt(2453750.5, 0.892862531);
        Assert.That(status, Is.EqualTo(0));
        Assert.That(tt1, Is.EqualTo(2453750.5).Within(1e-6));
        Assert.That(tt2, Is.EqualTo(0.8928551387488816828).Within(1e-12));
    }

    [Test]
    public void Tdbtcb_Tests() {
        var (status, tcb1, tcb2) = Sofa.Tdbtcb(2453750.5, 0.892855137);
        Assert.That(status, Is.EqualTo(0));
        Assert.That(tcb1, Is.EqualTo(2453750.5).Within(1e-6));
        Assert.That(tcb2, Is.EqualTo(0.8930195997253656716).Within(1e-12));
    }

    [Test]
    public void Tdbtt_Tests() {
        var (status, tt1, tt2) = Sofa.Tdbtt(2453750.5, 0.892855137, -0.000201);
        Assert.That(status, Is.EqualTo(0));
        Assert.That(tt1, Is.EqualTo(2453750.5).Within(1e-6));
        Assert.That(tt2, Is.EqualTo(0.8928551393263888889).Within(1e-12));
    }

    [Test]
    public void Tttai_Tests() {
        var (status, tai1, tai2) = Sofa.Tttai(2453750.5, 0.892482639);
        Assert.That(status, Is.EqualTo(0));
        Assert.That(tai1, Is.EqualTo(2453750.5).Within(1e-6));
        Assert.That(tai2, Is.EqualTo(0.892110139).Within(1e-12));
    }

    [Test]
    public void Tttcg_Tests() {
        var (status, tcg1, tcg2) = Sofa.Tttcg(2453750.5, 0.892482639);
        Assert.That(status, Is.EqualTo(0));
        Assert.That(tcg1, Is.EqualTo(2453750.5).Within(1e-6));
        Assert.That(tcg2, Is.EqualTo(0.8924900312508587113).Within(1e-12));
    }

    [Test]
    public void Tttdb_Tests() {
        var (status, tdb1, tdb2) = Sofa.Tttdb(2453750.5, 0.892855139, -0.000201);
        Assert.That(status, Is.EqualTo(0));
        Assert.That(tdb1, Is.EqualTo(2453750.5).Within(1e-6));
        Assert.That(tdb2, Is.EqualTo(0.8928551366736111111).Within(1e-12));
    }

    [Test]
    public void Ttut1_Tests() {
        var (status, ut11, ut12) = Sofa.Ttut1(2453750.5, 0.892855139, 64.8499);
        Assert.That(status, Is.EqualTo(0));
        Assert.That(ut11, Is.EqualTo(2453750.5).Within(1e-6));
        Assert.That(ut12, Is.EqualTo(0.8921045614537037037).Within(1e-12));
    }

    [Test]
    public void Ut1tai_Tests() {
        var (status, tai1, tai2) = Sofa.Ut1tai(2453750.5, 0.892104561, -32.6659);
        Assert.That(status, Is.EqualTo(0));
        Assert.That(tai1, Is.EqualTo(2453750.5).Within(1e-6));
        Assert.That(tai2, Is.EqualTo(0.8924826385462962963).Within(1e-12));
    }

    [Test]
    public void Ut1tt_Tests() {
        var (status, tt1, tt2) = Sofa.Ut1tt(2453750.5, 0.892104561, 64.8499);
        Assert.That(status, Is.EqualTo(0));
        Assert.That(tt1, Is.EqualTo(2453750.5).Within(1e-6));
        Assert.That(tt2, Is.EqualTo(0.8928551385462962963).Within(1e-12));
    }

    [Test]
    public void Ut1utc_Tests() {
        var (status, utc1, utc2) = Sofa.Ut1utc(2453750.5, 0.892104561, 0.3341);
        Assert.That(status, Is.EqualTo(0));
        Assert.That(utc1, Is.EqualTo(2453750.5).Within(1e-6));
        Assert.That(utc2, Is.EqualTo(0.8921006941018518519).Within(1e-12));
    }

    [Test]
    public void Utctai_Tests() {
        var (status, tai1, tai2) = Sofa.Utctai(2453750.5, 0.892100694);
        Assert.That(status, Is.EqualTo(0));
        Assert.That(tai1, Is.EqualTo(2453750.5).Within(1e-6));
        Assert.That(tai2, Is.EqualTo(0.8924826384444444444).Within(1e-12));
    }

    [Test]
    public void Utcut1_Tests() {
        var (status, ut11, ut12) = Sofa.Utcut1(2453750.5, 0.892100694, 0.3341);
        Assert.That(status, Is.EqualTo(0));
        Assert.That(ut11, Is.EqualTo(2453750.5).Within(1e-6));
        Assert.That(ut12, Is.EqualTo(0.8921045608981481481).Within(1e-12));
    }

    #endregion
    
    #region "Astronomy: Ecliptic Coordinates"
    /*######### Astronomy: Ecliptic Coordinates ##########################################################*/

    [Test]
    public void Eceq06_Tests() {
        double date1, date2, dl, db;
        date1 = 2456165.5;
        date2 = 0.401182685;
        dl = 5.1;
        db = -0.9;
        var (dr, dd) = Sofa.Eceq06(date1, date2, dl, db);
        Assert.That(dr, Is.EqualTo(5.533459733613627767).Within(1e-14));
        Assert.That(dd, Is.EqualTo(-1.246542932554480576).Within(1e-14));
    }

    [Test]
    public void Ecm06_Tests() {
        double date1, date2;
        date1 = 2456165.5;
        date2 = 0.401182685;
        var rm = Sofa.Ecm06(date1, date2);
        var rmExpected = 
            new double[,] {{0.9999952427708701137, -0.2829062057663042347e-2, -0.1229163741100017629e-2},
                           {0.3084546876908653562e-2, 0.9174891871550392514, 0.3977487611849338124},
                           {0.2488512951527405928e-5, -0.3977506604161195467, 0.9174935488232863071}};
        Assert.That(rm, Is.EqualTo(rmExpected).Within(1e-14));
    }

    [Test]
    public void Eqec06_Tests() {
        double date1, date2, dr, dd;
        date1 = 1234.5;
        date2 = 2440000.5;
        dr = 1.234;
        dd = 0.987;
        var (dl, db) = Sofa.Eqec06(date1, date2, dr, dd);
        Assert.That(dl, Is.EqualTo(1.342509918994654619).Within(1e-14));
        Assert.That(db, Is.EqualTo(0.5926215259704608132).Within(1e-14));
    }

    #endregion
    
    #region "Astronomy: Earth Rotation and Sidereal Time"
    /*######### Astronomy: Earth Rotation and Sidereal Time ##############################################*/

    [Test]
    public void Era00_Tests() {
        var era = Sofa.Era00(2400000.5, 54388.0);
        Assert.That(era, Is.EqualTo(0.4022837240028158102).Within(1e-12));
    }

    [Test]
    public void Gmst06_Tests() {
        var gmst = Sofa.Gmst06(2400000.5, 53736.0, 2400000.5, 53736.0);
        Assert.That(gmst, Is.EqualTo(1.754174971870091203).Within(1e-12));
    }

    [Test]
    public void Gst06a_Tests() {
        var gast = Sofa.Gst06a(2400000.5, 53736.0, 2400000.5, 53736.0);
        Assert.That(gast, Is.EqualTo(1.754166137675019159).Within(1e-12));
    }
    
    #endregion
    
    #region "Astronomy: Ephemerides"
    /*######### Astronomy: Ephemerides ###################################################################*/

    [Test]
    public void Epv00_Tests() {
        var (status, pvh, pvb) = Sofa.Epv00(2400000.5, 53411.52501161);
        Assert.That(status, Is.EqualTo(0));
        
        Assert.That(pvh.GetLength(0), Is.EqualTo(2));
        Assert.That(pvh.GetLength(1), Is.EqualTo(3));
        Assert.That(pvh[0,0], Is.EqualTo(-0.7757238809297706813).Within(1e-14));
        Assert.That(pvh[0,1], Is.EqualTo(0.5598052241363340596).Within(1e-14));
        Assert.That(pvh[0,2], Is.EqualTo(0.2426998466481686993).Within(1e-14));
        Assert.That(pvh[1,0], Is.EqualTo(-0.1091891824147313846e-1).Within(1e-15));
        Assert.That(pvh[1,1], Is.EqualTo(-0.1247187268440845008e-1).Within(1e-15));
        Assert.That(pvh[1,2], Is.EqualTo(-0.5407569418065039061e-2).Within(1e-15));
        
        Assert.That(pvb.GetLength(0), Is.EqualTo(2));
        Assert.That(pvb.GetLength(1), Is.EqualTo(3));
        Assert.That(pvb[0,0], Is.EqualTo(-0.7714104440491111971).Within(1e-14));
        Assert.That(pvb[0,1], Is.EqualTo(0.5598412061824171323).Within(1e-14));
        Assert.That(pvb[0,2], Is.EqualTo(0.2425996277722452400).Within(1e-14));
        Assert.That(pvb[1,0], Is.EqualTo(-0.1091874268116823295e-1).Within(1e-15));
        Assert.That(pvb[1,1], Is.EqualTo(-0.1246525461732861538e-1).Within(1e-15));
        Assert.That(pvb[1,2], Is.EqualTo(-0.5404773180966231279e-2).Within(1e-15));
    }

    [Test]
    public void Moon98_Tests() {
        var pv = Sofa.Moon98(2400000.5, 43999.9);
        var pvExpected = 
            new double[,] {{-0.2601295959971044180e-2, 0.6139750944302742189e-3, 0.2640794528229828909e-3},
                {-0.1244321506649895021e-3, -0.5219076942678119398e-3, -0.1716132214378462047e-3}};
        Assert.That(pv, Is.EqualTo(pvExpected).Within(1e-11));
    }

    [Test]
    public void Plan94_Tests() {
        // Test illegal planet number np=0:
        var (status, pv) = Sofa.Plan94(2400000.5, 1e6, 0);
        Assert.That(status, Is.EqualTo(-1));
        var pvExpected = new double[,] {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
        Assert.That(pv, Is.EqualTo(pvExpected));
        
        // Test illegal planet number np=10:
        (status, pv) = Sofa.Plan94(2400000.5, 1e6, 10);
        Assert.That(status, Is.EqualTo(-1));
        
        // Test year outside range:
        (status, pv) = Sofa.Plan94(2400000.5, -320000, 3);
        Assert.That(status, Is.EqualTo(1));
        pvExpected = new double[,] {{0.9308038666832975759, 0.3258319040261346000, 0.1422794544481140560}, 
            {-0.6429458958255170006e-2, 0.1468570657704237764e-1, 0.6406996426270981189e-2}};
        Assert.That(pv, Is.EqualTo(pvExpected).Within(1e-11));
        
        // Test normal case:
        (status, pv) = Sofa.Plan94(2400000.5, 43999.9, 1);
        Assert.That(status, Is.EqualTo(0));
        pvExpected = new double[,] {{0.2945293959257430832, -0.2452204176601049596, -0.1615427700571978153},
            {0.1413867871404614441e-1, 0.1946548301104706582e-1, 0.8929809783898904786e-2}};
        Assert.That(pv, Is.EqualTo(pvExpected).Within(1e-11));
    }

    #endregion
    
    
    #region "Astronomy: Fundamental Arguments"
    /*######### Astronomy: Fundamental Arguments #########################################################*/
    // No implemented functions to test.
    #endregion
    

    #region "Astronomy: Galactic Coordinates"
    /*######### Astronomy: Galactic Coordinates ##########################################################*/

    [Test]
    public void G2icrs_Tests() {
        double dl, db;
        dl =  5.5850536063818546461558105;
        db = -0.7853981633974483096156608;
        var (dr, dd) = Sofa.G2icrs(dl, db);
        Assert.That(dr, Is.EqualTo(5.9338074302227188048671).Within(1e-14));
        Assert.That(dd, Is.EqualTo(-1.1784870613579944551541).Within(1e-14));
    }

    [Test]
    public void Icrs2g_Tests() {
        double dr, dd;
        dr =  5.9338074302227188048671087;
        dd = -1.1784870613579944551540570;
        var (dl, db) = Sofa.Icrs2g(dr, dd);
        Assert.That(dl, Is.EqualTo(5.5850536063818546461558).Within(1e-14));
        Assert.That(db, Is.EqualTo(-0.7853981633974483096157).Within(1e-14));
    }
    
    #endregion
    
    #region "Astronomy: Geocentric/Geodetic Transformations"
    /*######### Astronomy: Geocentric/Geodetic Transformations ###########################################*/

    [Test]
    public void Gc2gd_Tests() {
        var xyz = new double[] {2e6, 3e6, 5.244e6};
        
        // Test illegal ellipsoid identifier n=0:
        var (status, elong, phi, height) = Sofa.Gc2gd(0, xyz);
        Assert.That(status, Is.EqualTo(-1));
        
        // Test normal case, WGS84:
        (status, elong, phi, height) = Sofa.Gc2gd(1, xyz); // n=1 signifies WGS84.
        Assert.That(status, Is.EqualTo(0));
        Assert.That(elong,  Is.EqualTo(0.9827937232473290680).Within(1e-14));
        Assert.That(phi,    Is.EqualTo(0.97160184819075459).Within(1e-14));
        Assert.That(height, Is.EqualTo(331.4172461426059892).Within(1e-8));

        // Test normal case, GRS80 ellipsoid:
        (status, elong, phi, height) = Sofa.Gc2gd(2, xyz); // n=2 signifies GRS80.
        Assert.That(status, Is.EqualTo(0));
        Assert.That(elong,  Is.EqualTo(0.9827937232473290680).Within(1e-14));
        Assert.That(phi,    Is.EqualTo(0.97160184820607853).Within(1e-14));
        Assert.That(height, Is.EqualTo(331.41731754844348).Within(1e-8));
        
        // Test normal case, WGS72 ellipsoid:
        (status, elong, phi, height) = Sofa.Gc2gd(3, xyz); // n=3 signifies WGS72.
        Assert.That(status, Is.EqualTo(0));
        Assert.That(elong,  Is.EqualTo(0.9827937232473290680).Within(1e-14));
        Assert.That(phi,    Is.EqualTo(0.9716018181101511937).Within(1e-14));
        Assert.That(height, Is.EqualTo(333.2770726130318123).Within(1e-8));
        
        // Test illegal ellipsoid identifier:
        (status, elong, phi, height) = Sofa.Gc2gd(4, xyz); // n not in range 1-3.
        Assert.That(status, Is.EqualTo(-1));
    }

    [Test]
    public void Gd2gc_Tests() {
        double elong = 3.1, phi = -0.5, height = 2500.0;
        
        // Test illegal ellipsoid identifier n=0:
        var (status, xyz) = Sofa.Gd2gc(0, elong, phi, height);
        Assert.That(status, Is.EqualTo(-1));
        
        // Test normal case, WGS84:
        (status, xyz) = Sofa.Gd2gc(1, elong, phi, height);
        Assert.That(status, Is.EqualTo(0));
        var xyzExpected = new double[] {-5599000.5577049947, 233011.67223479203, -3040909.4706983363};
        Assert.That(xyz,  Is.EqualTo(xyzExpected).Within(1e-7));
        
        // Test normal case, GRS80:
        (status, xyz) = Sofa.Gd2gc(2, elong, phi, height);
        Assert.That(status, Is.EqualTo(0));
        xyzExpected = new double[] {-5599000.5577260984, 233011.6722356702949, -3040909.4706095476};
        Assert.That(xyz,  Is.EqualTo(xyzExpected).Within(1e-7));
        
        // Test normal case, WGS72:
        (status, xyz) = Sofa.Gd2gc(3, elong, phi, height);
        Assert.That(status, Is.EqualTo(0));
        xyzExpected = new double[] {-5598998.7626301490, 233011.5975297822211, -3040908.6861467111};
        Assert.That(xyz,  Is.EqualTo(xyzExpected).Within(1e-7));
        
        // Test illegal ellipsoid identifier n=4:
        (status, xyz) = Sofa.Gd2gc(4, elong, phi, height);
        Assert.That(status, Is.EqualTo(-1));
    }

    #endregion
    
    #region "Astronomy: Gnomonic Projections"
    /*######### Astronomy: Gnomonic Projections ##########################################################*/

    [Test]
    public void Tpsts_Tests() {
        double xi, eta, raz, decz;
        xi = -0.03;
        eta = 0.07;
        raz = 2.3;
        decz = 1.5;
        var (ra, dec) = Sofa.Tpsts(xi, eta, raz, decz);
        Assert.That(ra,  Is.EqualTo(0.7596127167359629775).Within(1e-14));
        Assert.That(dec, Is.EqualTo(1.540864645109263028).Within(1e-13));
    }

    [Test]
    public void Tpstv_Tests() {
        double xi, eta, raz, decz;
        xi = -0.03;
        eta = 0.07;
        raz = 2.3;
        decz = 1.5;
        var vz = Sofa.S2c(raz, decz);
        var v = Sofa.Tpstv(xi, eta, vz);
        Assert.That(v.GetLength(0), Is.EqualTo(3));
        Assert.That(v[0],  Is.EqualTo(0.02170030454907376677).Within(1e-15));
        Assert.That(v[1],  Is.EqualTo(0.02060909590535367447).Within(1e-15));
        Assert.That(v[2],  Is.EqualTo(0.9995520806583523804).Within(1e-14)); 
    }

    [Test]
    public void Tpxes_Tests() {
        double ra, dec, raz, decz;
        ra = 1.3;
        dec = 1.55;
        raz = 2.3;
        decz = 1.5;
        var (status, xi, eta) = Sofa.Tpxes(ra, dec, raz, decz);
        Assert.That(status,  Is.EqualTo(0));
        Assert.That(xi,  Is.EqualTo(-0.01753200983236980595).Within(1e-15));
        Assert.That(eta, Is.EqualTo(0.05962940005778712891).Within(1e-15));
    }

    [Test]
    public void Tpxev_Tests() {
        double ra, dec, raz, decz;
        ra = 1.3;
        dec = 1.55;
        raz = 2.3;
        decz = 1.5;
        var v  = Sofa.S2c(ra, dec);
        var vz = Sofa.S2c(raz, decz);
        var (status, xi, eta) = Sofa.Tpxev(v, vz);
        Assert.That(status,  Is.EqualTo(0));
        Assert.That(xi,  Is.EqualTo(-0.01753200983236980595).Within(1e-15));
        Assert.That(eta, Is.EqualTo(0.05962940005778712891).Within(1e-15));
    }

    #endregion
    
    #region "Astronomy: Horizon/Equatorial Coordinates"
    /*######### Astronomy: Horizon/Equatorial Coordinates ################################################*/

    [Test]
    public void Ae2hd_Tests() {
        double az, el, phi;
        az = 5.5;
        el = 1.1;
        phi = 0.7;
        var (ha, dec) = Sofa.Ae2hd(az, el, phi);
        Assert.That(ha,  Is.EqualTo(0.5933291115507309663).Within(1e-14));
        Assert.That(dec, Is.EqualTo(0.9613934761647817620).Within(1e-14));
    }

    [Test]
    public void Hd2ae_Tests() {
        double ha, dec, phi;
        ha = 1.1;
        dec = 1.2;
        phi = 0.3;
        var (az, el) = Sofa.Hd2ae(ha, dec, phi);
        Assert.That(az, Is.EqualTo(5.916889243730066194).Within(1e-13));
        Assert.That(el, Is.EqualTo(0.4472186304990486228).Within(1e-14));
    }
    
    #endregion
    
    #region "Astronomy: Precession/Nutation/Polar Motion"
    /*######### Astronomy: Precession/Nutation/Polar Motion ##############################################*/

    [Test]
    public void S06a_Tests() {
        var s = Sofa.S06a(2400000.5, 52541.0);
        Assert.That(s, Is.EqualTo(-0.1340680437291812383e-7).Within(1e-18));
    }

    [Test]
    public void Xy06_Tests() {
        var (x, y) = Sofa.Xy06(2400000.5, 53736.0);
        Assert.That(x, Is.EqualTo(0.5791308486706010975e-3).Within(1e-15));
        Assert.That(y, Is.EqualTo(0.4020579816732958141e-4).Within(1e-16));
    }

    [Test]
    public void Xys06a_Tests() {
        var (x, y, s) = Sofa.Xys06a(2400000.5, 53736.0);
        Assert.That(x, Is.EqualTo(0.5791308482835292617e-3).Within(1e-14));
        Assert.That(y, Is.EqualTo(0.4020580099454020310e-4).Within(1e-15));
        Assert.That(s, Is.EqualTo(-0.1220032294164579896e-7).Within(1e-18));
    }
    
    #endregion
    
    #region "Astronomy: Star Catalog Conversions"
    /*######### Astronomy: Star Catalog Conversions ######################################################*/
    
    // No functions in this region were implemented.
    
    #endregion
    
    #region "Vector/Matrix: Initialization"
    /*######### Vector/Matrix: Initialization ############################################################*/
    
    // No functions in this region were implemented.
    
    #endregion
    
    #region "Vector/Matrix: Copy/Extend/Extract"
    /*######### Vector/Matrix: Copy/Extend/Extract #######################################################*/
        
    // No functions in this region were implemented.
    
    #endregion
    
    #region "Vector/Matrix: Build Rotations"
    /*######### Vector/Matrix: Build Rotations ###########################################################*/

    [Test]
    public void Rx_Tests() {
        double phi = 0.3456789;
        var r = new double[,] {{2, 3, 2}, {3, 2, 3}, {3, 4, 5}};
        var rotated = Sofa.Rx(phi, r);
        Assert.That(rotated[0,0], Is.EqualTo(2.0));
        Assert.That(rotated[0,1], Is.EqualTo(3.0));
        Assert.That(rotated[0,2], Is.EqualTo(2.0));
        Assert.That(rotated[1,0], Is.EqualTo(3.839043388235612460).Within(1e-12));
        Assert.That(rotated[1,1], Is.EqualTo(3.237033249594111899).Within(1e-12));
        Assert.That(rotated[1,2], Is.EqualTo(4.516714379005982719).Within(1e-12));
        Assert.That(rotated[2,0], Is.EqualTo(1.806030415924501684).Within(1e-12));
        Assert.That(rotated[2,1], Is.EqualTo(3.085711545336372503).Within(1e-12));
        Assert.That(rotated[2,2], Is.EqualTo(3.687721683977873065).Within(1e-12));
    }

    [Test]
    public void Ry_Tests() {
        double theta = 0.3456789;
        var r = new double[,] {{2, 3, 2}, {3, 2, 3}, {3, 4, 5}};
        var rotated = Sofa.Ry(theta, r);
        Assert.That(rotated[0,0], Is.EqualTo(0.8651847818978159930).Within(1e-12));
        Assert.That(rotated[0,1], Is.EqualTo(1.467194920539316554).Within(1e-12));
        Assert.That(rotated[0,2], Is.EqualTo(0.1875137911274457342).Within(1e-12));
        Assert.That(rotated[1,0], Is.EqualTo(3.0));
        Assert.That(rotated[1,1], Is.EqualTo(2.0));
        Assert.That(rotated[1,2], Is.EqualTo(3.0));
        Assert.That(rotated[2,0], Is.EqualTo(3.500207892850427330).Within(1e-12));
        Assert.That(rotated[2,1], Is.EqualTo(4.779889022262298150).Within(1e-12));
        Assert.That(rotated[2,2], Is.EqualTo(5.381899160903798712).Within(1e-12));
    }

    [Test]
    public void Rz_Tests() {
        double psi = 0.3456789;
        var r = new double[,] {{2, 3, 2}, {3, 2, 3}, {3, 4, 5}};
        var rotated = Sofa.Rz(psi, r);
        Assert.That(rotated[0,0], Is.EqualTo(2.898197754208926769).Within(1e-12));
        Assert.That(rotated[0,1], Is.EqualTo(3.500207892850427330).Within(1e-12));
        Assert.That(rotated[0,2], Is.EqualTo(2.898197754208926769).Within(1e-12));
        Assert.That(rotated[1,0], Is.EqualTo(2.144865911309686813).Within(1e-12));
        Assert.That(rotated[1,1], Is.EqualTo(0.865184781897815993).Within(1e-12));
        Assert.That(rotated[1,2], Is.EqualTo(2.144865911309686813).Within(1e-12));
        Assert.That(rotated[2,0], Is.EqualTo(3.0));
        Assert.That(rotated[2,1], Is.EqualTo(4.0));
        Assert.That(rotated[2,2], Is.EqualTo(5.0));
    }

    #endregion

    #region "Vector/Matrix: Spherical/Cartesian Conversions"
    /*######### Vector/Matrix: Spherical/Cartesian Conversions ###########################################*/

    [Test]
    public void S2c_Tests() {
        var c = Sofa.S2c(3.0123, -0.999);
        var cExpected = new double[] 
            {-0.5366267667260523906, 0.0697711109765145365, -0.8409302618566214041};
        Assert.That(c, Is.EqualTo(cExpected).Within(1e-12));
    }

    [Test]
    public void C2s_Tests() {
        var p = new double[] {100, -50, 25};
        var (theta, phi) = Sofa.C2s(p);
        Assert.That(theta, Is.EqualTo(-0.4636476090008061162).Within(1e-14));
        Assert.That(phi,   Is.EqualTo(0.2199879773954594463).Within(1e-14));
    }

    [Test]
    public void S2p_Tests() {
        var p = Sofa.S2p(-3.21, 0.123, 0.456);
        var pExpected = new double[] 
            {-0.4514964673880165228, 0.0309339427734258688, 0.0559466810510877933};
        Assert.That(p, Is.EqualTo(pExpected).Within(1e-12));
    }

    [Test]
    public void P2s_Tests() {
        var p = new double[] {100, -50, 25};
        var (theta, phi, r) = Sofa.P2s(p);
        Assert.That(theta, Is.EqualTo(-0.4636476090008061162).Within(1e-12));
        Assert.That(phi,   Is.EqualTo(0.2199879773954594463).Within(1e-12));
        Assert.That(r,     Is.EqualTo(114.5643923738960002).Within(1e-9));
    }

    [Test]
    public void S2pv_Tests() {
        var pv = Sofa.S2pv(-3.21, 0.123, 0.456, -7.8e-6, 9.01e-6, -1.23e-5);
        Assert.That(pv.GetLength(0), Is.EqualTo(2));
        Assert.That(pv.GetLength(1), Is.EqualTo(3));
        Assert.That(pv[0,0], Is.EqualTo(-0.4514964673880165228).Within(1e-12));
        Assert.That(pv[0,1], Is.EqualTo(0.0309339427734258688).Within(1e-12));
        Assert.That(pv[0,2], Is.EqualTo(0.0559466810510877933).Within(1e-12));
        Assert.That(pv[1,0], Is.EqualTo(0.1292270850663260170e-4).Within(1e-16));
        Assert.That(pv[1,1], Is.EqualTo(0.2652814182060691422e-5).Within(1e-16));
        Assert.That(pv[1,2], Is.EqualTo(0.2568431853930292259e-5).Within(1e-16));
    }

    [Test]
    public void Pv2s_Tests() {
        var pv = new double[,] {{-0.4514964673880165, 0.03093394277342585, 0.05594668105108779},
            {1.292270850663260e-5, 2.652814182060692e-6, 2.568431853930293e-6}};
        var (theta, phi, r, td, pd, rd) = Sofa.Pv2s(pv);
        Assert.That(theta, Is.EqualTo(3.073185307179586515).Within(1e-12));
        Assert.That(phi,   Is.EqualTo(0.1229999999999999992).Within(1e-12));
        Assert.That(r,     Is.EqualTo(0.4559999999999999757).Within(1e-12));
        Assert.That(td,    Is.EqualTo(-0.7800000000000000364e-5).Within(1e-16));
        Assert.That(pd,    Is.EqualTo(0.9010000000000001639e-5).Within(1e-16));
        Assert.That(rd,    Is.EqualTo(-0.1229999999999999832e-4).Within(1e-16));
    }

    #endregion
    
    #region "Vector/Matrix: Operations on Vectors"
    /*######### Vector/Matrix: Operations on Vectors #####################################################*/
    
    // No functions in this region were implemented.
    
    #endregion
    
    #region "Vector/Matrix: Operations on Matrices"
    /*######### Vector/Matrix: Operations on Matrices ####################################################*/

    // TODO: Resume here...
    
    [Test]
    public void Rxr_Tests() { }

    [Test]
    public void Tr_Tests() { }
    
    
    #endregion
    
    
    
}