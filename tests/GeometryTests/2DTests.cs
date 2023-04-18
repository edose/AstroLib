using System.Diagnostics.CodeAnalysis;
using NUnit.Framework;
using AstroLib.Core.Geometry;

// ReSharper disable CompareOfFloatsByEqualityOperator
// ReSharper disable CommentTypo
// ReSharper disable InconsistentNaming
// ReSharper disable IdentifierTypo
// ReSharper disable InlineTemporaryVariable
// ReSharper disable JoinDeclarationAndInitializer

namespace AstroLibTests.GeometryTests;

[TestFixture]
[SuppressMessage("Assertion", "NUnit2045:Use Assert.Multiple")]
public class Point2DTests {
    
    private readonly Point2D p1 = new Point2D(3, 4);
    private readonly Point2D p2 = new Point2D(6, 8);
    private readonly Point2D p3 = new Point2D(11, 11);
    private readonly Vector2D v1 = new Vector2D(2.5, 3.7);    
    
    // [SetUp]
    // public void Setup() { }
    
    [Test]
    public void Constructor_Tests() {
        // // Default constructor:
        // var pt1 = new Point2D {X = 33, Y = 44};
        // Assert.That(pt1.X, Is.EqualTo(33.0));
        // Assert.That(pt1.Y, Is.EqualTo(44.0));
        
        // From 2 doubles:
        var pt2 = new Point2D(34, 45);
        Assert.That(pt2.X, Is.EqualTo(34.0));
        Assert.That(pt2.Y, Is.EqualTo(45.0));
        
        // From Tuple of 2 doubles:
        var tuple3 = Tuple.Create(35.0, 46.0); 
        var pt3 = new Point2D(tuple3);
        Assert.That(Tuple.Create(pt3.X, pt3.Y), Is.EqualTo(tuple3));
        
        // From Tuple of 2 ints:
        var tuple3i = Tuple.Create(42, 53); 
        var pt3i = new Point2D(tuple3i);
        Assert.That(pt3i.X, Is.EqualTo(42.0));
        Assert.That(pt3i.Y, Is.EqualTo(53.0));
        
        // From array of 2 doubles:
        var array4 = new double[] {36, 47};
        var pt4 = new Point2D(array4);
        Assert.That(pt4.X, Is.EqualTo(36.0));
        Assert.That(pt4.Y, Is.EqualTo(47.0));
    }

    [Test]
    public void Method_Tests() {
        Point2D pResult;
        Vector2D vResult;
        // .Add(vector):
        pResult = p1.Add(v1);
        Assert.That(pResult.X, Is.EqualTo(p1.X + v1.Dx));
        Assert.That(pResult.Y, Is.EqualTo(p1.Y + v1.Dy));
        
        // .VectorTo(point):
        vResult = p1.VectorTo(p3);
        Assert.That(vResult.Dx, Is.EqualTo(p3.X - p1.X));
        Assert.That(vResult.Dy, Is.EqualTo(p3.Y - p1.Y));
        
        // .VectorFrom(point):
        vResult = p1.VectorFrom(p3);
        Assert.That(vResult.Dx, Is.EqualTo(p1.X - p3.X));
        Assert.That(vResult.Dy, Is.EqualTo(p1.Y - p3.Y));
        
        // Test consistency on reversal:
        Assert.That(p1.VectorTo(p2), Is.EqualTo(p2.VectorFrom(p1)));
        Assert.That(p1.VectorTo(p2), Is.EqualTo(p2.VectorTo(p1).Reversed));
        
        // .DistanceToLine(point, point) where points define the line:
        var linePointA = new Point2D(0, 3);
        var linePointB = new Point2D(4, 0);
        var point = new Point2D(4, 3);
        var distance = point.DistanceToLine(linePointA, linePointB);
        Assert.That(distance, Is.EqualTo(2.4).Within(1E-9));
        
        var distanceLinePointsReversed = point.DistanceToLine(linePointB, linePointA);
        Assert.That(distance, Is.EqualTo(distanceLinePointsReversed).Within(1E-9));
    }
}

[TestFixture]
[SuppressMessage("Assertion", "NUnit2045:Use Assert.Multiple")]
public class Vector2DTests {

    private readonly Vector2D v1 = new Vector2D(2.5, 3.7);
    private readonly Vector2D v2 = new Vector2D(3.5, 4.7);
    private readonly Vector2D v3 = new Vector2D(6.2, 8.9);
    private readonly Point2D p1 = new Point2D(3, 4);

    // [SetUp]
    // public void Setup() { }

    [Test]
    public void Constructor_Tests() {
        // // Default constructor:
        // var vector1 = new Vector2D {Dx = 33, Dy = 44};
        // Assert.That(vector1.Dx, Is.EqualTo(33.0));
        // Assert.That(vector1.Dy, Is.EqualTo(44.0));
        
        // From doubles:
        var vector2 = new Vector2D(34, 22);
        Assert.That(vector2.Dx, Is.EqualTo(34.0));
        Assert.That(vector2.Dy, Is.EqualTo(22.0));
        
        // From Tuple of 2 doubles:
        var tuple = Tuple.Create(32.0, 43.0);
        var vector3 = new Vector2D(tuple);
        Assert.That(vector3.Dx, Is.EqualTo(32.0));
        Assert.That(vector3.Dy, Is.EqualTo(43.0));
        
        // From Tuple of 2 ints:
        var tupleInt = Tuple.Create(52, 63);
        var vector3Int = new Vector2D(tupleInt);
        Assert.That(vector3Int.Dx, Is.EqualTo(52.0));
        Assert.That(vector3Int.Dy, Is.EqualTo(63.0));
        
        // From array of 2 doubles:
        var array = new double[] {21, 32};
        var vector4 = new Vector2D(array);
        Assert.That(vector4.Dx, Is.EqualTo(21.0));
        Assert.That(vector4.Dy, Is.EqualTo(32.0));
    }

    [Test]
    public void Property_Tests() {
        Assert.That(v1.Dx, Is.EqualTo(2.5));
        Assert.That(v1.Dy, Is.EqualTo(3.7));
        Assert.That(v1.Length2, Is.EqualTo(2.5*2.5 + 3.7*3.7).Within(1E-9));
        Assert.That(v1.Length, Is.EqualTo(Math.Sqrt(v1.Length2)).Within(1E-9));
        Assert.That(v1.Direction, Is.EqualTo(Math.Atan2(3.7, 2.5)).Within(1E-9));
        var v1Reversed = v1.Reversed;
        Assert.That(v1Reversed.Dx, Is.EqualTo(-2.5));
        Assert.That(v1Reversed.Dy, Is.EqualTo(-3.7));
    }
    
    [Test]
    public void Method_Tests() {
        Vector2D vResult;
        
        // .Subtract(vector):
        vResult = v1.Subtract(v3);
        Assert.That(vResult.Dx, Is.EqualTo(-3.7).Within(1E-9));
        Assert.That(vResult.Dy, Is.EqualTo(-5.2).Within(1E-9));
        
        // .MultiplyBy(factor):
        vResult = v2.MultiplyBy(5.0);
        Assert.That(vResult.Dx, Is.EqualTo(5.0 * v2.Dx).Within(1E-9));
        Assert.That(vResult.Dy, Is.EqualTo(5.0 * v2.Dy).Within(1E-9));
        
        // .DivideBy(factor):
        vResult = v2.DivideBy(4.5);
        Assert.That(vResult.Dx, Is.EqualTo(v2.Dx / 4.5).Within(1E-9));
        Assert.That(vResult.Dy, Is.EqualTo(v2.Dy / 4.5).Within(1E-9));
        
        // .DotProduct(vector):
        var dotResult = v2.DotProduct(v3);
        Assert.That(dotResult, Is.EqualTo(v2.Dx * v3.Dx + v2.Dy * v3.Dy).Within(1E-9));
        
        // .AngleWith(vector):
        var angleResult = v2.AngleWith(v3);
        var angleExpected = v3.Direction - v2.Direction;
        Assert.That(angleResult, Is.EqualTo(angleExpected).Within(1E-9));
        Assert.That(v2.AngleWith(v2), Is.EqualTo(0.0));
        Assert.That(v2.AngleWith(v2.Reversed), Is.EqualTo(Math.PI).Within(1E-9));
    }
}

[TestFixture]
[SuppressMessage("Assertion", "NUnit2045:Use Assert.Multiple")]
public class Rectangle2DTests {

    private readonly Point2D a = new Point2D(0, 3);
    private readonly Point2D b = new Point2D(4, 0);
    private readonly Point2D c = new Point2D(10, 8);

    // [SetUp]
    // public void Setup() { }

    [Test]
    public void Constructor_Tests() {
        // Default constructor:
        var r = new Rectangle2D(a, b, c);
        Assert.That(r.IsValid, Is.EqualTo(true));
        Assert.That(r.A.X, Is.EqualTo(0.0));
        Assert.That(r.C.Y, Is.EqualTo(8.0));
        Assert.That(r.D.X, Is.EqualTo(6.0));
        Assert.That(r.D.Y, Is.EqualTo(11.0));
        }

    [Test]
    public void Property_Tests() {
        var r = new Rectangle2D(a, b, c);
        Assert.That(r.Area, Is.EqualTo(50.0).Within(1E-9));
        Assert.That(r.IsValid, Is.EqualTo(true));
    }

    [Test]
    public void Method_Tests() {
        var r = new Rectangle2D(a, b, c);
        
        // .ContainsPoint(point):
        var pointInside = new Point2D(1.0, 3.0);
        Assert.That(r.ContainsPoint(pointInside, includeEdges: true), Is.EqualTo(true));
        Assert.That(r.ContainsPoint(pointInside, includeEdges: false), Is.EqualTo(true));
        var pointCorner = new Point2D(0.0, 3.0);
        Assert.That(r.ContainsPoint(pointCorner, includeEdges: true), Is.EqualTo(true));
        Assert.That(r.ContainsPoint(pointCorner, includeEdges: false), Is.EqualTo(false));
        var pointEdge = new Point2D(3.0, 7.0);
        Assert.That(r.ContainsPoint(pointEdge, includeEdges: true), Is.EqualTo(true));
        Assert.That(r.ContainsPoint(pointEdge, includeEdges: false), Is.EqualTo(false));
        var pointOutside = new Point2D(20.0, 30.0);
        Assert.That(r.ContainsPoint(pointOutside, includeEdges: true), Is.EqualTo(false));
        Assert.That(r.ContainsPoint(pointOutside, includeEdges: false), Is.EqualTo(false));
    }
}

[TestFixture]
[SuppressMessage("Assertion", "NUnit2045:Use Assert.Multiple")]
public class Circle2DTests {

    private readonly Point2D origin = new Point2D(10, 8);

    // [SetUp]
    // public void Setup() { }

    [Test]
    public void Constructor_Tests() {
        // Default constructor:
        var circle = new Circle2D(origin, 5.0);
        Assert.That(circle.Origin, Is.EqualTo(origin));
        Assert.That(circle.Radius, Is.EqualTo(5.0));
    }
    
    [Test]
    public void Property_Tests() {
        var circle = new Circle2D(origin, 25.0);
        Assert.That(circle.Origin, Is.EqualTo(origin));
        Assert.That(circle.Radius, Is.EqualTo(25.0));
        Assert.That(circle.Area, Is.EqualTo(Math.PI * 25.0 * 25.0).Within(1E-9));
    }
    
    [Test]
    public void Method_Tests() {
        var circle = new Circle2D(origin, 15.0);
        
        // .ContainsPoint(point):
        var pointInside = new Point2D(11.0, 9.0);
        Assert.That(circle.ContainsPoint(pointInside, includeEdges: true), Is.EqualTo(true));
        Assert.That(circle.ContainsPoint(pointInside, includeEdges: false), Is.EqualTo(true));
        var pointEdge = new Point2D(-2.0, -1.0);
        Assert.That(circle.ContainsPoint(pointEdge, includeEdges: true), Is.EqualTo(true));
        Assert.That(circle.ContainsPoint(pointEdge, includeEdges: false), Is.EqualTo(false));
        var pointOutside = new Point2D(-3.0, -2.0);
        Assert.That(circle.ContainsPoint(pointOutside, includeEdges: true), Is.EqualTo(false));
        Assert.That(circle.ContainsPoint(pointOutside, includeEdges: false), Is.EqualTo(false));
    }
}
