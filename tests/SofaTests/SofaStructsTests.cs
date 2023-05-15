using AstroLib.Core;
using NUnit.Framework;
using AstroLib.Core.SOFA;
// ReSharper disable IdentifierTypo

namespace AstroLibTests.SofaTests; 

[TestFixture]
public class SofaStructsTests {

    [Test]
    public void iauASTROM_Tests() {
        var iau = new Sofa.iauASTROM {
            pmt = 10.0,
            eb = new double[3] {20, 21, 22},
            eh = new double[3] {30, 31, 32},
            em = 41,
            v = new double[3] {50, 51, 52},
            bm1 = 61,
            bpn = new double[9] {70, 71, 72, 73, 74, 75, 76, 77, 78},
            along = 81,
            phi = 87,
            xpl = 91,
            ypl = 101,
            sphi = 111,
            cphi = 121,
            diurab = 131,
            eral = 141,
            refa = 151,
            refb = 161
        };
        Assert.That(iau.pmt, Is.EqualTo(10));
        Assert.That(iau.eb, Is.EqualTo(new double[3] {20, 21, 22}));
        Assert.That(iau.eh, Is.EqualTo(new double[3] {30, 31, 32}));
        Assert.That(iau.em, Is.EqualTo(41));
        Assert.That(iau.v, Is.EqualTo(new double[3] {50, 51, 52}));
        Assert.That(iau.bm1, Is.EqualTo(61));
        Assert.That(iau.bpn, Is.EqualTo(new double[9] {70, 71, 72, 73, 74, 75, 76, 77, 78}));
        Assert.That(iau.along, Is.EqualTo(81));
        Assert.That(iau.phi, Is.EqualTo(87));
        Assert.That(iau.xpl, Is.EqualTo(91));
        Assert.That(iau.ypl, Is.EqualTo(101));
        Assert.That(iau.sphi, Is.EqualTo(111));
        Assert.That(iau.cphi, Is.EqualTo(121));
        Assert.That(iau.diurab, Is.EqualTo(131));
        Assert.That(iau.eral, Is.EqualTo(141));
        Assert.That(iau.refa, Is.EqualTo(151));
        Assert.That(iau.refb, Is.EqualTo(161));
    }

    [Test]
    public void TranslateToCSharp_Tests() {
        var cpp = new Sofa.iauASTROM {
            pmt = 10.0,
            eb = new double[3] {20, 21, 22},
            eh = new double[3] {30, 31, 32},
            em = 41,
            v = new double[3] {50, 51, 52},
            bm1 = 61,
            bpn = new double[9] {70, 71, 72, 73, 74, 75, 76, 77, 78},
            along = 81,
            phi = 87,
            xpl = 91,
            ypl = 101,
            sphi = 111,
            cphi = 121,
            diurab = 131,
            eral = 141,
            refa = 151,
            refb = 161
        };
        var csharp = Sofa.TranslateToCsharp(cpp);
        Assert.That(csharp.pmt, Is.EqualTo(10));
        Assert.That(csharp.eb, Is.EqualTo(new double[3] {20, 21, 22}));
        Assert.That(csharp.eh, Is.EqualTo(new double[3] {30, 31, 32}));
        Assert.That(csharp.em, Is.EqualTo(41));
        Assert.That(csharp.v, Is.EqualTo(new double[3] {50, 51, 52}));
        Assert.That(csharp.bm1, Is.EqualTo(61));
        Assert.That(csharp.bpn, Is.EqualTo(new double[3,3] {{70, 71, 72}, {73, 74, 75}, {76, 77, 78}}));
        Assert.That(csharp.along, Is.EqualTo(81));
        Assert.That(csharp.phi, Is.EqualTo(87));
        Assert.That(csharp.xpl, Is.EqualTo(91));
        Assert.That(csharp.ypl, Is.EqualTo(101));
        Assert.That(csharp.sphi, Is.EqualTo(111));
        Assert.That(csharp.cphi, Is.EqualTo(121));
        Assert.That(csharp.diurab, Is.EqualTo(131));
        Assert.That(csharp.eral, Is.EqualTo(141));
        Assert.That(csharp.refa, Is.EqualTo(151));
        Assert.That(csharp.refb, Is.EqualTo(161));
    }

}