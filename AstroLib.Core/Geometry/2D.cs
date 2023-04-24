// ReSharper disable MemberCanBePrivate.Global

using System.Diagnostics.CodeAnalysis;

namespace AstroLib.Core.Geometry;

/// <summary>Represents one (X,Y) point on a plane. X and Y are dimensionless.</summary>
public class Point2D : IEquatable<Point2D> {
    
    public double X { get; init; }
    public double Y { get; init; }

    /// <summary>Private default constructor, to prevent user construction without input data.</summary>
    private Point2D() {}

    /// <summary>Constructor from 2 doubles, x and y.</summary>
    public Point2D(double x, double y) { X = x; Y = y; }

    /// <summary>Constructor from Tuple of 2 doubles.</summary>
    public Point2D(Tuple<double, double>xy) { (X, Y) = xy; }
    
    /// <summary>Constructor from Tuple of 2 ints.</summary>
    public Point2D(Tuple<int, int>xy) { (X, Y) = xy; }

    /// <summary>Constructor from array.</summary>
    public Point2D(double[] xy) { X = xy[0]; Y = xy[1]; }
    

    /// <summary> Returns true iff other Point2D equals this one, exactly.</summary> 
    [SuppressMessage("ReSharper", "CompareOfFloatsByEqualityOperator")]
    public bool Equals(Point2D? other) {
        if (other == null) { return false; }
        return (other.X == X) && (other.Y == Y);
    }

    /// <summary>Return a new point displaced from this point by a vector.</summary>
    public Point2D Add(Vector2D dxy) {
        return new Point2D {X = X + dxy.Dx, Y = Y + dxy.Dy};
    }

    /// <summary>Conceptual alias for .Add(). </summary>
    public Point2D DisplaceBy(Vector2D dxy) { return this.Add(dxy); }

    /// <summary> Return a new vector from this point to other point. </summary>
    public Vector2D VectorTo(Point2D other) {
        return new Vector2D(this, other);
    }

    /// <summary> Return a new vector from other point to this point. </summary>
    public Vector2D VectorFrom(Point2D other) {
        return new Vector2D(other, this);
    }

    /// <summary>Return distance from this point to another point.</summary>
    public double DistanceTo(Point2D other) {
        return Math.Sqrt((X - other.X) * (X - other.X) + (Y - other.Y) * (Y - other.Y));
    }

    /// <summary> Return distance from this point to a line defined by two points.</summary>
    public double DistanceToLine(Point2D a, Point2D b) {
        var distAb = new Vector2D(a, b).Length;
        if (distAb == 0) {
            return (new Vector2D(this, a)).Length;
        }
        return Math.Abs((b.X - a.X) * (a.Y - Y) - (b.Y - a.Y) * (a.X - X)) / distAb;
    }

    /// <summary> Return distance from this point to a line *segment* defined by two points.</summary>
    public double DistanceToSegment(Point2D a, Point2D b) {
        var ax = new Vector2D(a, this);
        var ab = new Vector2D(a, b);
        var dpAxAb = ax.DotProduct(ab);
        var bx = new Vector2D(b, this);
        var ba = new Vector2D(b, a);
        var dpBxBa = bx.DotProduct(ba);
        
        // If point a and point b are directly on, or are on the the same side of the normal
        // from this point to the line:
        if (Math.Sign(dpAxAb) == Math.Sign(dpBxBa)) {
            return Math.Min(ax.Length, bx.Length);
        }        
        // Point a and point must be on opposite sides of the normal from this point to the line:
        return this.DistanceToLine(a, b);
    }
}

/// <summary>Represents a vector in a plane. Dx and Dy are dimensionless.</summary>
public class Vector2D : IEquatable<Vector2D> {
    public double Dx { get; init; }
    public double Dy { get; init; }
    public double Length2 => this.Dx * this.Dx + this.Dy * this.Dy;
    public double Length => Math.Sqrt(this.Length2);
    public double Direction => Math.Atan2(this.Dy, this.Dx);
    public Vector2D Reversed => new Vector2D(-Dx, -Dy);

    /// <summary>Private default constructor, to prevent user construction without input data.</summary>
    private Vector2D() {}

    /// <summary>Constructor from 2 doubles, dx and dy.</summary>
    public Vector2D(double dx, double dy) { Dx = dx; Dy = dy; }
    
    /// <summary>Constructor from two Point2D objects.</summary>
    public Vector2D(Point2D from, Point2D to) { Dx = to.X - from.X; Dy = to.Y - from.Y; }
    
    /// <summary>Constructor from Tuple of 2 doubles.</summary>
    public Vector2D(Tuple<double, double>dxy) { (Dx, Dy) = dxy; }
    
    /// <summary>Constructor from Tuple of 2 ints.</summary>
    public Vector2D(Tuple<int, int>dxy) { (Dx, Dy) = dxy; }
    
    /// <summary>Constructor from array.</summary>
    public Vector2D(double[] dxy) { Dx = dxy[0]; Dy = dxy[1]; }

    [SuppressMessage("ReSharper", "CompareOfFloatsByEqualityOperator")]
    public bool Equals(Vector2D? other) {
        if (other == null) { return false; }
        return (other.Dx == Dx) && (other.Dy == Dy);
    }

    /// <summary>Add two vectors, returning the sum vector.</summary>
    public Vector2D Add(Vector2D other) {
        return new Vector2D {Dx = Dx + other.Dx, Dy = Dy + other.Dy};
    }
    
    /// <summary>Subtract other vector from this vector, returning the difference vector.</summary>
    public Vector2D Subtract(Vector2D other) {
        return new Vector2D {Dx = Dx - other.Dx, Dy = Dy - other.Dy};
    }

    /// <summary>Multiply this vector by a scalar factor.</summary>
    public Vector2D MultiplyBy(double factor) {
        return new Vector2D {Dx = factor * Dx, Dy = factor * Dy};
    }
    
    /// <summary>Divide this vector by a scalar factor.</summary>
    public Vector2D DivideBy(double divisor) {
        return new Vector2D {Dx = Dx / divisor, Dy = Dy / divisor};
    }
    
    /// <summary>Return dot product between this vector and another.</summary>
    public double DotProduct(Vector2D other) {
        return Dx * other.Dx + Dy * other.Dy;
    }
    
    /// <summary>Return angle (radians) between this vector and another vector.</summary>
    public double AngleWith(Vector2D other) {
        return Math.Acos(DotProduct(other) / (Length * other.Length));
    }
}

/// <summary>Represents a rectangle on a plane.</summary>
public class Rectangle2D {
    public Point2D A { get; init; }
    public Point2D B { get; init; }
    public Point2D C { get; init; }
    public Point2D D { get; init; }
    private readonly Vector2D ab, bc;
    public readonly double Area;
    public readonly bool IsValid;

    /// <summary>Main constructor. The three points must be any three consecutive corners of the rectangle,
    /// i.e., of a right triangle forming half the rectangle. Points may be passed in
    /// clockwise or counterclockwise order.</summary>
    public Rectangle2D(Point2D a, Point2D b, Point2D c) {
        ab = a.VectorTo(b);
        bc = b.VectorTo(c);
        A = a; B = b; C = c;
        D = c.Add(ab.Reversed);
        Area = ab.Length * bc.Length;
        IsValid = (ab.Length > 0) &&
                  (bc.Length > 0) &&
                  (Math.Abs(ab.DotProduct(bc)) < 1.0E-9 * Math.Min(ab.Length, bc.Length));
    }

    /// <summary>Returns true if this rectangle contains this Point2D, else return false.</summary>
    /// <param name="xy">The point to be tested (Point2D object).</param>
    /// <param name="includeEdges">True if a point directly on the rectangle's boundaries is to be
    /// considered contained, False if it is to be considered not contained.</param>
    /// <returns>True iff this rectangle contains this point.</returns>
    public bool ContainsPoint(Point2D xy, bool includeEdges=true) {
        var dotAbAb = ab.Length2;
        var dotBcBc = bc.Length2;
        var dotAbPt = ab.DotProduct(A.VectorTo(xy));
        var dotBcPt = bc.DotProduct(B.VectorTo(xy));
        if (includeEdges)
            return (0 <= dotAbPt) && (dotAbPt <= dotAbAb) && (0 <= dotBcPt) && (dotBcPt <= dotBcBc);
        return (0 < dotAbPt) && (dotAbPt < dotAbAb) && (0 < dotBcPt) && (dotBcPt < dotBcBc);
    }
}

/// <summary>Represents a circle on Cartesian plane.</summary>
public class Circle2D {
    public Point2D Origin { get; init; }
    public double Radius { get; init; }
    public double Area { get; init; }

    /// <summary>Constructor, given origin point and radius.</summary>
    public Circle2D(Point2D origin, double radius) {
        Origin = origin;
        Radius = radius;
        Area = Math.PI * radius * radius;
    }

    /// <summary>Returns true if this circle contains this Point2D, else return false.</summary>
    /// <param name="xy">The point to be tested (Point2D object).</param>
    /// <param name="includeEdges">True if a point directly on the circle's boundary is to be
    /// considered contained, False if it is to be considered not contained.</param>
    /// <returns>True iff this circle contains this point.</returns>
    public bool ContainsPoint(Point2D xy, bool includeEdges = true) {
        var distance2 = Origin.VectorTo(xy).Length2;
        // Console.WriteLine($"distance2:{distance2}  " +
        //                   $"radius^2:{Radius * Radius}  " +
        //                   $"diff:{distance2 - Radius * Radius}");
        if (includeEdges)
            return (distance2 <= Radius * Radius);
        return (distance2 < Radius * Radius);
    }
}
