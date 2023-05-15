using System.Diagnostics.CodeAnalysis;
using AstroLib.Core;
using NUnit.Framework;


// ReSharper disable CompareOfFloatsByEqualityOperator
// ReSharper disable CommentTypo
// ReSharper disable InconsistentNaming
// ReSharper disable IdentifierTypo
// ReSharper disable InlineTemporaryVariable
// ReSharper disable JoinDeclarationAndInitializer

namespace AstroLibTests.CoreTests;

[TestFixture]
[SuppressMessage("Assertion", "NUnit2045:Use Assert.Multiple")]
public class AstroMathTests {
    
    // [SetUp]
    // public void Setup() { }

    [Test]
    public void EuclidianModuloTests() {
        Assert.That(AstroMath.EuclidianModulo(7.5, 4), Is.EqualTo(3.5));
        Assert.That(AstroMath.EuclidianModulo(7.5, -4), Is.EqualTo(3.5));
        Assert.That(AstroMath.EuclidianModulo(-7.5, 4), Is.EqualTo(0.5));
        Assert.That(AstroMath.EuclidianModulo(-7.5, -4), Is.EqualTo(0.5));
        Assert.That(AstroMath.EuclidianModulo(0, 4), Is.EqualTo(0.0));
        Assert.That(AstroMath.EuclidianModulo(0, -4), Is.EqualTo(0.0));
        Assert.That(AstroMath.EuclidianModulo(103, 4), Is.EqualTo(3));
        Assert.That(AstroMath.EuclidianModulo(103, -4), Is.EqualTo(3));
        Assert.That(AstroMath.EuclidianModulo(7.5, 0), Is.Null);
        Assert.That(AstroMath.EuclidianModulo(-7.5, 0), Is.Null);
        Assert.That(AstroMath.EuclidianModulo(0.3333, 1.777), Is.EqualTo(0.3333)); // exactness check.
    }

    [Test]
    public void WrapTests() {
        Assert.That(AstroMath.Wrap(-455, 0, 360), Is.EqualTo(-455 - (-720)));
        Assert.That(AstroMath.Wrap(-55, 0, 360), Is.EqualTo(360 - 55));
        Assert.That(AstroMath.Wrap(0, 0, 360), Is.EqualTo(0));
        Assert.That(AstroMath.Wrap(132, 0, 360), Is.EqualTo(132));
        Assert.That(AstroMath.Wrap(359.999, 0, 360), Is.EqualTo(359.999));
        Assert.That(AstroMath.Wrap(360, 0, 360), Is.EqualTo(0));
        Assert.That(AstroMath.Wrap(367, 0, 360), Is.EqualTo(7));
        Assert.That(AstroMath.Wrap(667, 0, 360), Is.EqualTo(667 - 360));
        Assert.That(AstroMath.Wrap(1.3333, 0.342, 1.866), 
            Is.EqualTo(1.3333));  // exactness check.
        
        Assert.That(AstroMath.Wrap(-455, -180, 180), Is.EqualTo(-455 - (-360)));
        Assert.That(AstroMath.Wrap(-255, -180, 180), Is.EqualTo(-255 - (-360)));
        Assert.That(AstroMath.Wrap(-180.125, -180, 180), Is.EqualTo(-180.125 + 360));
        Assert.That(AstroMath.Wrap(-180, -180, 180), Is.EqualTo(-180));
        Assert.That(AstroMath.Wrap(-55, -180, 180), Is.EqualTo(-55));
        Assert.That(AstroMath.Wrap(0, -180, 180), Is.EqualTo(0));
        Assert.That(AstroMath.Wrap(-57, -180, 180), Is.EqualTo(-57));
        Assert.That(AstroMath.Wrap(179.875, -180, 180), Is.EqualTo(179.875));
        Assert.That(AstroMath.Wrap(180, -180, 180), Is.EqualTo(-180));
        Assert.That(AstroMath.Wrap(180.125, -180, 180), Is.EqualTo(180.125 - 360));
        Assert.That(AstroMath.Wrap(567, -180, 180), Is.EqualTo(567 - 720));
    }

    [Test]
    public void ClampTests() {
        Assert.That(AstroMath.Clamp(-111, -23, 31), Is.EqualTo(-23));
        Assert.That(AstroMath.Clamp(-23, -23, 31), Is.EqualTo(-23));
        Assert.That(AstroMath.Clamp(0, -23, 31), Is.EqualTo(0));
        Assert.That(AstroMath.Clamp(31, -23, 31), Is.EqualTo(31));
        Assert.That(AstroMath.Clamp(110, -23, 31), Is.EqualTo(31));
        Assert.That(AstroMath.Clamp(1.3333, 0.342, 1.866), 
            Is.EqualTo(1.3333)); // exactness check
        
        Assert.That(AstroMath.Clamp(-111, null, 31), Is.EqualTo(-111));
        Assert.That(AstroMath.Clamp(-23, null, 31), Is.EqualTo(-23));
        Assert.That(AstroMath.Clamp(0, null, 31), Is.EqualTo(0));
        Assert.That(AstroMath.Clamp(31, null, 31), Is.EqualTo(31));
        Assert.That(AstroMath.Clamp(110, null, 31), Is.EqualTo(31));
        
        Assert.That(AstroMath.Clamp(-111, -23, null), Is.EqualTo(-23));
        Assert.That(AstroMath.Clamp(-23, -23, null), Is.EqualTo(-23));
        Assert.That(AstroMath.Clamp(0, -23, null), Is.EqualTo(0));
        Assert.That(AstroMath.Clamp(31, -23, null), Is.EqualTo(31));
        Assert.That(AstroMath.Clamp(110, -23, null), Is.EqualTo(110));
        
        Assert.That(AstroMath.Clamp(-111, null, null), Is.EqualTo(-111));
        Assert.That(AstroMath.Clamp(-23, null, null), Is.EqualTo(-23));
        Assert.That(AstroMath.Clamp(0, null, null), Is.EqualTo(0));
        Assert.That(AstroMath.Clamp(31, null, null), Is.EqualTo(31));
        Assert.That(AstroMath.Clamp(110, null, null), Is.EqualTo(110));
    }
    
    
    [Test]
    public void IsInRangeTests() {
        Assert.That(AstroMath.IsInRange(-111, -23, 31), Is.False);
        Assert.That(AstroMath.IsInRange(-23, -23, 31), Is.True);
        Assert.That(AstroMath.IsInRange(0, -23, 31), Is.True);
        Assert.That(AstroMath.IsInRange(31, -23, 31), Is.True);
        Assert.That(AstroMath.IsInRange(110, -23, 31), Is.False);
        
        Assert.That(AstroMath.IsInRange(-111, null, 31), Is.True);
        Assert.That(AstroMath.IsInRange(-23, null, 31), Is.True);
        Assert.That(AstroMath.IsInRange(0, null, 31), Is.True);
        Assert.That(AstroMath.IsInRange(31, null, 31), Is.True);
        Assert.That(AstroMath.IsInRange(110.0, null, 31), Is.False);
        
        Assert.That(AstroMath.IsInRange(-111, -23, null), Is.False);
        Assert.That(AstroMath.IsInRange(-23, -23, null), Is.True);
        Assert.That(AstroMath.IsInRange(0, -23, null), Is.True);
        Assert.That(AstroMath.IsInRange(31, -23, null), Is.True);
        Assert.That(AstroMath.IsInRange(110, -23, null), Is.True);
        
        Assert.That(AstroMath.IsInRange(-111, null, null), Is.True);
        Assert.That(AstroMath.IsInRange(-23, null, null), Is.True);
        Assert.That(AstroMath.IsInRange(0, null, null), Is.True);
        Assert.That(AstroMath.IsInRange(31, null, null), Is.True);
        Assert.That(AstroMath.IsInRange(110, null, null), Is.True);
    }

    [Test]
    public void Render1dTo2dArrayTests() {
        var array1d = new double[6] {3, 2, 6, 5, 9, 8};
        var array2d = AstroMath.Reshape1dTo2dArray(array1d, 2, 3);
        var array2dExpected = new double[,] {{3, 2, 6}, {5, 9, 8} };
        Assert.That(array2d, Is.EqualTo(array2dExpected));
        
        // Test throwing exception for mismatched sizes:
        Assert.That(() => AstroMath.Reshape1dTo2dArray(array1d, 2, 2), 
            Throws.TypeOf<ArgumentException>());
    }
}