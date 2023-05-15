using NUnit.Framework;
using AstroLib.Core.SOFA;
// ReSharper disable IdentifierTypo
// ReSharper disable RedundantExplicitArrayCreation

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

    #region "Astronomy: Astrometry"
    /*######### Astronomy: Astrometry ####################################################################*/

    [Test]
    public void Ab_Tests() {
        // Adapted to AstroLib from http://iausofa.org/2021_0512_C/sofa/t_sofa_c.html.
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
        var bpnExpected = new double[,] {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};
        Assert.That(astrom.bpn, Is.EqualTo(bpnExpected));
        
        // TODO: Resume here...

        
    }
    
    public void Apci13_Tests(){ }
    
    public void Apco13_Tests(){ }
    
    public void Apcs13_Tests() {}
    
    public void Aper13_Tests() {}
    
    public void Apio13_Tests() {}
    
    public void Atcc13_Tests() {}
    
    public void Atccq_Tests() {}
    
    public void Atci13_Tests() {}
    
    public void Atciq_Tests() {}
    
    public void Atciqz_Tests() {}
    
    public void Atco13_Tests() {}
    
    public void Atic13_Tests() {}
    
    public void Aticq_Tests() {}
    
    public void Atio13_Tests() {}
    
    public void Atioq_Tests() {}
    
    public void Atoc13_Tests() {}
    
    public void Atoi13_Tests() {}
    
    public void Atoiq_Tests() {}
    
    public void Pmsafe_Tests() {}
    
    public void Pvtob_Tests() {}
    
    public void Refco_Tests() {}
    
    public void Starpm_Tests() {}

    #endregion

}