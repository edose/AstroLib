using System.Diagnostics.CodeAnalysis;
using AstroLib.Core.Geometry;
using NUnit.Framework;

// ReSharper disable CompareOfFloatsByEqualityOperator
// ReSharper disable CommentTypo
// ReSharper disable InconsistentNaming
// ReSharper disable IdentifierTypo
// ReSharper disable InlineTemporaryVariable
// ReSharper disable JoinDeclarationAndInitializer

namespace AstroLibTests.GeometryTests;

[TestFixture]
[SuppressMessage("Assertion", "NUnit2045:Use Assert.Multiple")]
public class Point3DTests {
    
    private readonly Point3D p1 = new Point3D(3, 4, 6);
    private readonly Point3D p2 = new Point3D(6, 8, 11);
    private readonly Point3D p3 = new Point3D(11, 11, 32);
    private readonly Vector3D v1 = new Vector3D(2.5, 3.7, 5.1);
    
    // [SetUp]
    // public void Setup() { }

    [Test]
    public void Constructor_Tests() {
        // Default constructor:
        var pt1 = new Point3D {X = 33, Y = 44, Z = 56};
        Assert.That(pt1.X, Is.EqualTo(33.0));
        Assert.That(pt1.Y, Is.EqualTo(44.0));
        Assert.That(pt1.Z, Is.EqualTo(56.0));
        
        // From doubles:
        var pt2 = new Point3D(34, 45, 73);
        Assert.That(pt2.X, Is.EqualTo(34.0));
        Assert.That(pt2.Y, Is.EqualTo(45.0));
        Assert.That(pt2.Z, Is.EqualTo(73.0));
        
        // From Tuple of 3 doubles:
        var tuple3 = Tuple.Create(35.0, 46.0, 21.0); 
        var pt3 = new Point3D(tuple3);
        Assert.That(Tuple.Create(pt3.X, pt3.Y, pt3.Z), Is.EqualTo(tuple3));

        // From array of 3 doubles:
        var array4 = new double[] {36, 47, 32};
        var pt4 = new Point3D(array4);
        Assert.That(pt4.X, Is.EqualTo(36.0));
        Assert.That(pt4.Y, Is.EqualTo(47.0));
        Assert.That(pt4.Z, Is.EqualTo(32.0));
    }

    [Test]
    public void Method_Tests() {
        Point3D pResult;
        Vector3D vResult;
        
        // .Add(vector):
        pResult = p1.Add(v1);
        Assert.That(pResult.X, Is.EqualTo(p1.X + v1.Dx));
        Assert.That(pResult.Y, Is.EqualTo(p1.Y + v1.Dy));
        Assert.That(pResult.Z, Is.EqualTo(p1.Z + v1.Dz));
        
        // .VectorTo(point):
        vResult = p1.VectorTo(p3);
        Assert.That(vResult.Dx, Is.EqualTo(p3.X - p1.X));
        Assert.That(vResult.Dy, Is.EqualTo(p3.Y - p1.Y));
        Assert.That(vResult.Dz, Is.EqualTo(p3.Z - p1.Z));

        // .VectorFrom(point):
        vResult = p1.VectorFrom(p3);
        Assert.That(vResult.Dx, Is.EqualTo(p1.X - p3.X));
        Assert.That(vResult.Dy, Is.EqualTo(p1.Y - p3.Y));
        Assert.That(vResult.Dz, Is.EqualTo(p1.Z - p3.Z));

        // Test consistency on reversal:
        Assert.That(p1.VectorTo(p2), Is.EqualTo(p2.VectorFrom(p1)));
    }
}

[TestFixture]
[SuppressMessage("Assertion", "NUnit2045:Use Assert.Multiple")]
public class Vector3DTests {

    private readonly Vector3D v1 = new Vector3D(2.5, 3.7, 4.1);
    private readonly Vector3D v2 = new Vector3D(3.5, 4.7, 6.6);
    private readonly Vector3D v3 = new Vector3D(6.2, 8.9, 5.9);
    private readonly Point3D p1 = new Point3D(3, 4, 8);

    // [SetUp]
    // public void Setup() { }

    [Test]
    public void Constructor_Tests() {
        // Default constructor:
        var vector1 = new Vector3D {Dx = 33, Dy = 44, Dz = 51};
        Assert.That(vector1.Dx, Is.EqualTo(33.0));
        Assert.That(vector1.Dy, Is.EqualTo(44.0));
        Assert.That(vector1.Dz, Is.EqualTo(51.0));

        // From doubles:
        var vector2 = new Vector3D(34, 22, 13);
        Assert.That(vector2.Dx, Is.EqualTo(34.0));
        Assert.That(vector2.Dy, Is.EqualTo(22.0));
        Assert.That(vector2.Dz, Is.EqualTo(13.0));
        
        // From Tuple of 3 doubles:
        var tuple = Tuple.Create(32.0, 43.0, 53.0);
        var vector3 = new Vector3D(tuple);
        Assert.That(vector3.Dx, Is.EqualTo(32.0));
        Assert.That(vector3.Dy, Is.EqualTo(43.0));
        Assert.That(vector3.Dz, Is.EqualTo(53.0));
        
        // From array of 3 doubles:
        var array = new double[] {21, 32, 72};
        var vector4 = new Vector3D(array);
        Assert.That(vector4.Dx, Is.EqualTo(21.0));
        Assert.That(vector4.Dy, Is.EqualTo(32.0));
        Assert.That(vector4.Dz, Is.EqualTo(72.0));
        }

    [Test]
    public void Property_Tests() {
        Assert.That(v1.Dx, Is.EqualTo(2.5));
        Assert.That(v1.Dy, Is.EqualTo(3.7));
        Assert.That(v1.Dz, Is.EqualTo(4.1));
        Assert.That(v1.Length2, Is.EqualTo(2.5*2.5 + 3.7*3.7 + 4.1 * 4.1).Within(1E-9));
        Assert.That(v1.Length, Is.EqualTo(Math.Sqrt(v1.Length2)).Within(1E-9));
        var v1Reversed = v1.Reversed;
        Assert.That(v1Reversed.Dx, Is.EqualTo(-2.5));
        Assert.That(v1Reversed.Dy, Is.EqualTo(-3.7));
        Assert.That(v1Reversed.Dz, Is.EqualTo(-4.1));
    }

    [Test]
    public void Method_Tests() {
        Vector3D vResult;

        // .Add(vector):
        vResult = v2.Add(v1);
        Assert.That(vResult.Dx, Is.EqualTo(6.0).Within(1E-9));
        Assert.That(vResult.Dy, Is.EqualTo(8.4).Within(1E-9));
        Assert.That(vResult.Dz, Is.EqualTo(10.7).Within(1E-9));
        
        // .Subtract(vector):
        vResult = v3.Subtract(v1);
        Assert.That(vResult.Dx, Is.EqualTo(3.7).Within(1E-9));
        Assert.That(vResult.Dy, Is.EqualTo(5.2).Within(1E-9));
        Assert.That(vResult.Dz, Is.EqualTo(1.8).Within(1E-9));
        
        // .MultiplyBy(factor):
        vResult = v2.MultiplyBy(5.0);
        Assert.That(vResult.Dx, Is.EqualTo(5.0 * v2.Dx).Within(1E-9));
        Assert.That(vResult.Dy, Is.EqualTo(5.0 * v2.Dy).Within(1E-9));
        Assert.That(vResult.Dz, Is.EqualTo(5.0 * v2.Dz).Within(1E-9));
        
        // .DivideBy(divisor):
        vResult = v2.DivideBy(5.0);
        Assert.That(vResult.Dx, Is.EqualTo(v2.Dx / 5.0).Within(1E-9));
        Assert.That(vResult.Dy, Is.EqualTo(v2.Dy / 5.0).Within(1E-9));
        Assert.That(vResult.Dz, Is.EqualTo(v2.Dz / 5.0).Within(1E-9));
        
        // .DotProduct(vector):
        var dotResult = v2.DotProduct(v3);
        Assert.That(dotResult, Is.EqualTo(v2.Dx * v3.Dx + v2.Dy * v3.Dy + v2.Dz * v3.Dz).Within(1E-9));
        
        // .AngleWith(vector):
        Assert.That(v1.AngleWith(v1), Is.EqualTo(0.0));
        Assert.That(v1.AngleWith(v1.Reversed), Is.EqualTo(Math.PI).Within(1E-9));
        var va = new Vector3D(1, 1, 0);
        var vb = new Vector3D(-1, 1, -1);
        Assert.That(va.AngleWith(vb), Is.EqualTo(Math.PI / 2.0).Within(1E-9));
        Assert.That(va.AngleWith(vb), Is.EqualTo(vb.AngleWith(va)));
    }


}
