// ReSharper disable IdentifierTypo
// ReSharper disable MemberCanBePrivate.Global

namespace AstroLib.Core.Geometry; 

/// <summary>Represents a point in three-dimensional space.</summary>
public record Point3D {
    
    public double X { get; init; }
    public double Y { get; init; }
    public double Z { get; init; }

    /// <summary>Default constructor.</summary>
    public Point3D() {}

    /// <summary>Constructor from 3 doubles.</summary>
    public Point3D(double x, double y, double z) { X = x; Y = y; Z = z; }

    /// <summary>Constructor from Tuple of 3 doubles.</summary>
    public Point3D(Tuple<double, double, double>xyz) { (X, Y, Z) = xyz; }

    /// <summary>Constructor from array of 3 doubles.</summary>
    public Point3D(double[] xyz) { X = xyz[0]; Y = xyz[1]; Z = xyz[2]; }
    
    /// <summary>Return a new point displaced from this point by a vector.</summary>
    public Point3D Add(Vector3D dxyz) {
        return new Point3D {X = this.X + dxyz.Dx, Y = this.Y + dxyz.Dy, Z = this.Z + dxyz.Dz};
    }
    
    // /// <summary>Subtract another point from this one, returning a vector from other to this.</summary>
    // public Vector3D Subtract(Point3D other) {
    //     return new Vector3D(other, this);
    // }   
    //
    // /// <summary>Subtract a vector from this point, returning a new point displaced by
    // /// the vector's inverse. </summary>
    // public Point3D Subtract(Vector3D dxyz){
    //     return new Point3D {X = this.X - dxyz.Dx, Y = this.Y - dxyz.Dy, Z = this.Z - dxyz.Dz};
    // }
    
    /// <summary> Return a new vector from this point to other point. </summary>
    public Vector3D VectorTo(Point3D other) {
        return new Vector3D(this, other);
    }
    
    /// <summary> Return a new vector from other point to this point. </summary>
    public Vector3D VectorFrom(Point3D other) {
        return new Vector3D(other, this);
    }
}

/// <summary>Represents a vector in three-dimensional space.</summary>
public record Vector3D {
    public double Dx { get; init; }
    public double Dy { get; init; }
    public double Dz { get; init; }
    public double Length2 => this.Dx * this.Dx + this.Dy * this.Dy + this.Dz * this.Dz;
    public double Length => Math.Sqrt(this.Length2);
    public Vector3D Reversed => new Vector3D(-Dx, -Dy, -Dz);

    /// <summary>Default constructor.</summary>
    public Vector3D() {}

    /// <summary>Constructor from 3 doubles.</summary>
    public Vector3D(double dx, double dy, double dz) { Dx = dx; Dy = dy; Dz = dz; }

    /// <summary>Constructor from two Point3D objects.</summary>
    public Vector3D(Point3D from, Point3D to) {
        Dx = to.X - from.X;
        Dy = to.Y - from.Y;
        Dz = to.Z - from.Z;
    }

    /// <summary>Constructor from Tuple.</summary>
    public Vector3D(Tuple<double, double, double>dxyz) { (Dx, Dy, Dz) = dxyz; }
    
    /// <summary>Constructor from array.</summary>
    public Vector3D(double[] dxyz) { Dx = dxyz[0]; Dy = dxyz[1]; Dz = dxyz[2]; }
    
    /// <summary>Add two vectors, returning the sum vector.</summary>
    public Vector3D Add(Vector3D other) {
        return new Vector3D {Dx = this.Dx + other.Dx, Dy = this.Dy + other.Dy, Dz = this.Dz + other.Dz};
    }
    
    /// <summary>Subtract other vector from this vector, returning the difference vector.</summary>
    public Vector3D Subtract(Vector3D other) {
        return new Vector3D {Dx = this.Dx - other.Dx, Dy = this.Dy - other.Dy, Dz = this.Dz - other.Dz};
    }
    
    /// <summary>Multiply this vector by a scalar factor.</summary>
    public Vector3D MultiplyBy(double factor) {
        return new Vector3D {Dx = factor * this.Dx, Dy = factor * this.Dy, Dz = factor * this.Dz};
    }
    
    /// <summary>Divide this vector by a scalar factor.</summary>
    public Vector3D DivideBy(double divisor) {
        return new Vector3D {Dx = this.Dx / divisor, Dy = this.Dy / divisor, Dz = this.Dz / divisor};
    }
    
    /// <summary>Return dot product between this vector and another.</summary>
    public double DotProduct(Vector3D other) {
        return this.Dx * other.Dx + this.Dy * other.Dy + this.Dz * other.Dz;
    }
    
    /// <summary>Return angle (radians) between this vector and another vector.</summary>
    public double AngleWith(Vector3D other) {
        return Math.Acos(this.DotProduct(other) / (this.Length * other.Length));
    }   
}
    
