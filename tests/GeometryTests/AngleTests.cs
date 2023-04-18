using NUnit.Framework;
using AstroLib.Core;
using AstroLib.Core.Geometry;

namespace AstroLibTests.GeometryTests;

[TestFixture]
public class AngleClassTests {

    [Test]
    public void Constructor_Tests() {
        // From radians:
        var a1 = new Angle(0.375);
        Assert.IsInstanceOf(typeof(Angle), a1);
        Assert.That(a1.Radians, Is.EqualTo(0.375));

        // From degrees:
        var a2 = Angle.FromDegrees(38.8);
        Assert.IsInstanceOf(typeof(Angle), a2);
        Assert.That(a2.Radians, Is.EqualTo(38.8 / Constants.DegreesPerRadian).Within(1E-9));
    }

    [Test]
    public void Property_Tests() {
        var a1 = new Angle(0.375);
        Assert.That(a1.Radians, Is.EqualTo(0.375));
        Assert.That(a1.Degrees, Is.EqualTo(0.375 * Constants.DegreesPerRadian).Within(1E-9));
        Assert.That(a1.Arcminutes, Is.EqualTo(60.0 * 0.375 * Constants.DegreesPerRadian).Within(1E-9));
        Assert.That(a1.Arcseconds, Is.EqualTo(3600.0 * 0.375 * Constants.DegreesPerRadian).Within(1E-9));
        // Assert.That(a1.InputMin, Is.Null);
        // Assert.That(a1.InputMax, Is.Null);
    }

    [Test]
    public void Method_Tests() {
        // .Equals():
        var a1 = new Angle(0.375);
        var a2 = new Angle(0.375);
        var a3 = new Angle(0.111);
        Assert.That(a1, Is.EqualTo(a2));
        Assert.That(a1, Is.Not.EqualTo(a3));
        
        // .ToString():
        var words = a1.ToString().Split();
        Assert.That(words[0], Is.EqualTo("Angle:"));
        Assert.That(words[1], Is.EqualTo("0.375"));
        Assert.That(words[2], Is.EqualTo("radians"));
        Assert.That(words[3], Is.EqualTo("="));
        Assert.That(words[4].StartsWith("21.4859"), Is.True);
        Assert.That(words[5].StartsWith("degrees"));
    }
}

[TestFixture]
public class LatitudeClassTests {

    [Test]
    public void Constructor_Tests() {
        // From radians:
        var lat1 = new Latitude(0.375);
        Assert.IsInstanceOf(typeof(Latitude), lat1);
        Assert.That(lat1.Radians, Is.EqualTo(0.375));

        // From degrees:
        var lat2 = Latitude.FromDegrees(38.8);
        Assert.IsInstanceOf(typeof(Latitude), lat2);
        Assert.That(lat2.Radians, Is.EqualTo(38.8 / Constants.DegreesPerRadian).Within(1E-9));
        
        // From input outside limit of [-90, +90] degrees:
        Assert.That(() => Latitude.FromDegrees(-91), Throws.TypeOf<ArgumentException>());
        Assert.That(() => Latitude.FromDegrees(+91), Throws.TypeOf<ArgumentException>());
        Assert.That(Latitude.FromDegrees(-91, true), Is.EqualTo(Latitude.FromDegrees(-90)));
        Assert.That(Latitude.FromDegrees(+91, true), Is.EqualTo(Latitude.FromDegrees(+90)));
    }
    
    [Test]
    public void Property_Tests() {
        var lat1 = new Latitude(0.375);
        Assert.That(lat1.Radians, Is.EqualTo(0.375));
        Assert.That(lat1.Degrees, Is.EqualTo(0.375 * Constants.DegreesPerRadian).Within(1E-9));
        Assert.That(lat1.Arcminutes, Is.EqualTo(60.0 * 0.375 * Constants.DegreesPerRadian).Within(1E-9));
        Assert.That(lat1.Arcseconds, Is.EqualTo(3600.0 * 0.375 * Constants.DegreesPerRadian).Within(1E-9));
        Assert.That(lat1.ClampMin, Is.EqualTo(-Math.PI / 2.0));
        Assert.That(lat1.ClampMax, Is.EqualTo(+Math.PI / 2.0));
    }
    
    [Test]
    public void Method_Tests() {
        // .Equals():
        var lat1 = new Latitude(0.375);
        var lat2 = new Latitude(0.375);
        var lat3 = new Latitude(0.111);
        Assert.That(lat1, Is.EqualTo(lat2));
        Assert.That(lat1, Is.Not.EqualTo(lat3));
        
        // .ToPolarAngle():
        var pa1 = new Latitude(0.375).ToPolarAngle();
        Assert.IsInstanceOf(typeof(PolarAngle), pa1);
        Assert.That(pa1.Radians, Is.EqualTo(Math.PI / 2.0 - 0.375).Within(1E-9));

        // .ToString():
        var words = lat1.ToString().Split();
        Assert.That(words[0], Is.EqualTo("Latitude:"));
        Assert.That(words[1], Is.EqualTo("0.375"));
        Assert.That(words[2], Is.EqualTo("radians"));
        Assert.That(words[3], Is.EqualTo("="));
        Assert.That(words[4].StartsWith("21.4859"), Is.True);
        Assert.That(words[5].StartsWith("degrees"));
    }
}

[TestFixture]
public class PolarAngleClassTests {

    [Test]
    public void Constructor_Tests() {
        // From radians:
        var pa1 = new PolarAngle(0.375);
        Assert.IsInstanceOf(typeof(PolarAngle), pa1);
        Assert.That(pa1.Radians, Is.EqualTo(0.375));

        // From degrees:
        var pa2 = PolarAngle.FromDegrees(38.8);
        Assert.IsInstanceOf(typeof(PolarAngle), pa2);
        Assert.That(pa2.Radians, Is.EqualTo(38.8 / Constants.DegreesPerRadian).Within(1E-9));
        
        // From input outside limit of [0, +180] degrees:
        Assert.That(() => PolarAngle.FromDegrees(-11), Throws.TypeOf<ArgumentException>());
        Assert.That(() => PolarAngle.FromDegrees(+191), Throws.TypeOf<ArgumentException>());
        Assert.That(PolarAngle.FromDegrees(-11, true), 
            Is.EqualTo(PolarAngle.FromDegrees(0)));
        Assert.That(PolarAngle.FromDegrees(+191, true), 
            Is.EqualTo(PolarAngle.FromDegrees(+180, true)));
    }

    [Test]
    public void Property_Tests() {
        var pa1 = new PolarAngle(0.375);
        Assert.That(pa1.Radians, Is.EqualTo(0.375));
        Assert.That(pa1.Degrees, Is.EqualTo(0.375 * Constants.DegreesPerRadian).Within(1E-9));
        Assert.That(pa1.Arcminutes, Is.EqualTo(60.0 * 0.375 * Constants.DegreesPerRadian).Within(1E-9));
        Assert.That(pa1.Arcseconds, Is.EqualTo(3600.0 * 0.375 * Constants.DegreesPerRadian).Within(1E-9));
        Assert.That(pa1.ClampMin, Is.EqualTo(0));
        Assert.That(pa1.ClampMax, Is.EqualTo(Math.PI));
    }
    
    [Test]
    public void Method_Tests() {
        // .Equals():
        var pa1 = new PolarAngle(0.375);
        var pa2 = new PolarAngle(0.375);
        var pa3 = new PolarAngle(0.111);
        Assert.That(pa1, Is.EqualTo(pa2));
        Assert.That(pa1, Is.Not.EqualTo(pa3));
        
        // .ToLatitude():
        var lat1 = new PolarAngle(0.375).ToLatitude();
        Assert.IsInstanceOf(typeof(Latitude), lat1);
        Assert.That(lat1.Radians, Is.EqualTo(Math.PI / 2.0 - 0.375).Within(1E-9));

        // .ToString():
        var words = pa1.ToString().Split();
        Assert.That(words[0], Is.EqualTo("PolarAngle:"));
        Assert.That(words[1], Is.EqualTo("0.375"));
        Assert.That(words[2], Is.EqualTo("radians"));
        Assert.That(words[3], Is.EqualTo("="));
        Assert.That(words[4].StartsWith("21.4859"), Is.True);
        Assert.That(words[5].StartsWith("degrees"));
    }
}

[TestFixture]
public class LongitudeClassTests {

    [Test]
    public void Constructor_Tests() {
        // From radians:
        var lon1 = new Longitude(0.375);
        Assert.IsInstanceOf(typeof(Longitude), lon1);
        Assert.That(lon1.Radians, Is.EqualTo(0.375));

        // From degrees:
        var lon2 = Longitude.FromDegrees(38.8);
        Assert.IsInstanceOf(typeof(Longitude), lon2);
        Assert.That(lon2.Radians, Is.EqualTo(38.8 / Constants.DegreesPerRadian).Within(1E-9));

        // From input outside limit of [-180, +180] degrees:
        Assert.That(() => Longitude.FromDegrees(-184), Throws.TypeOf<ArgumentException>());
        Assert.That(() => Longitude.FromDegrees(+185), Throws.TypeOf<ArgumentException>());
        Assert.That(Longitude.FromDegrees(-184, true).Radians,
            Is.EqualTo(Longitude.FromDegrees(-184 + 360).Radians).Within(1E-9));
        Assert.That(Longitude.FromDegrees(+191, true).Radians,
            Is.EqualTo(Longitude.FromDegrees(+191 - 360, true).Radians).Within(1E-8));
    }
    
    [Test]
    public void Property_Tests() {
        var lon1 = new Longitude(0.875);
        Assert.That(lon1.Radians, Is.EqualTo(0.875));
        Assert.That(lon1.Degrees, Is.EqualTo(0.875 * Constants.DegreesPerRadian).Within(1E-9));
        Assert.That(lon1.Arcminutes, Is.EqualTo(60.0 * 0.875 * Constants.DegreesPerRadian).Within(1E-9));
        Assert.That(lon1.Arcseconds, Is.EqualTo(3600.0 * 0.875 * Constants.DegreesPerRadian).Within(1E-9));
        Assert.That(lon1.WrapMin, Is.EqualTo(-Math.PI));
        Assert.That(lon1.WrapMax, Is.EqualTo(+Math.PI));
    }
    
    [Test]
    public void Method_Tests() {
        // .Equals():
        var lon1 = new Longitude(0.875);
        var lon2 = new Longitude(0.875);
        var lon3 = new Longitude(0.111);
        Assert.That(lon1, Is.EqualTo(lon2));
        Assert.That(lon1, Is.Not.EqualTo(lon3));

        // .ToString():
        var words = lon1.ToString().Split();
        Assert.That(words[0], Is.EqualTo("Longitude:"));
        Assert.That(words[1], Is.EqualTo("0.875"));
        Assert.That(words[2], Is.EqualTo("radians"));
        Assert.That(words[3], Is.EqualTo("="));
        Assert.That(words[4].StartsWith("50.1338"), Is.True);
        Assert.That(words[5].StartsWith("degrees"));
    }
}

[TestFixture]
public class DeclinationClassTests {

    [Test]
    public void Constructor_Tests() {
        // From radians:
        var dec1 = new Declination(0.375);
        Assert.IsInstanceOf(typeof(Declination), dec1);
        Assert.That(dec1.Radians, Is.EqualTo(0.375));

        // From degrees:
        var dec2 = Declination.FromDegrees(38.8);
        Assert.IsInstanceOf(typeof(Declination), dec2);
        Assert.That(dec2.Radians, Is.EqualTo(38.8 / Constants.DegreesPerRadian).Within(1E-9));
        
                
        // From input outside limit of [-90, +90] degrees:
        Assert.That(() => Declination.FromDegrees(-91), Throws.TypeOf<ArgumentException>());
        Assert.That(() => Declination.FromDegrees(+91), Throws.TypeOf<ArgumentException>());
        Assert.That(Declination.FromDegrees(-91, true), Is.EqualTo(Declination.FromDegrees(-90)));
        Assert.That(Declination.FromDegrees(+91, true), Is.EqualTo(Declination.FromDegrees(+90)));
    }
    
    [Test]
    public void Property_Tests() {
        var lat1 = new Declination(0.375);
        Assert.That(lat1.Radians, Is.EqualTo(0.375));
        Assert.That(lat1.Degrees, Is.EqualTo(0.375 * Constants.DegreesPerRadian).Within(1E-9));
        Assert.That(lat1.Arcminutes, Is.EqualTo(60.0 * 0.375 * Constants.DegreesPerRadian).Within(1E-9));
        Assert.That(lat1.Arcseconds, Is.EqualTo(3600.0 * 0.375 * Constants.DegreesPerRadian).Within(1E-9));
        Assert.That(lat1.WrapMin, Is.EqualTo(-Math.PI / 2.0));
        Assert.That(lat1.WrapMax, Is.EqualTo(+Math.PI / 2.0));
    }
    
    [Test]
    public void Method_Tests() {
        // .Equals():
        var dec1 = new Declination(0.375);
        var dec2 = new Declination(0.375);
        var dec3 = new Declination(0.111);
        Assert.That(dec1, Is.EqualTo(dec2));
        Assert.That(dec1, Is.Not.EqualTo(dec3));
        
        // .ToPolarAngle():
        var pa1 = new Declination(0.375).ToPolarAngle();
        Assert.IsInstanceOf(typeof(PolarAngle), pa1);
        Assert.That(pa1.Radians, Is.EqualTo(Math.PI / 2.0 - 0.375).Within(1E-9));

        // .ToString():
        var words = dec1.ToString().Split();
        Assert.That(words[0], Is.EqualTo("Declination:"));
        Assert.That(words[1], Is.EqualTo("0.375"));
        Assert.That(words[2], Is.EqualTo("radians"));
        Assert.That(words[3], Is.EqualTo("="));
        Assert.That(words[4].StartsWith("21.4859"), Is.True);
        Assert.That(words[5].StartsWith("degrees"));
    }
}

[TestFixture]
public class RightAscensionClassTests {

    [Test]
    public void Constructor_Tests() {
        // From radians:
        var ra1 = new RightAscension(0.375);
        Assert.IsInstanceOf(typeof(RightAscension), ra1);
        Assert.That(ra1.Radians, Is.EqualTo(0.375));

        // From degrees:
        var ra2 = RightAscension.FromDegrees(38.8);
        Assert.IsInstanceOf(typeof(RightAscension), ra2);
        Assert.That(ra2.Radians, Is.EqualTo(38.8 / Constants.DegreesPerRadian).Within(1E-9));

        // From input outside range of [0, 360) degrees:
        Assert.That(() => RightAscension.FromDegrees(-14), Throws.TypeOf<ArgumentException>());
        Assert.That(() => RightAscension.FromDegrees(+385), Throws.TypeOf<ArgumentException>());
        Assert.That(RightAscension.FromDegrees(-14, true).Radians,
            Is.EqualTo(RightAscension.FromDegrees(-14 + 360).Radians).Within(1E-9));
        Assert.That(RightAscension.FromDegrees(+385, true).Radians,
            Is.EqualTo(RightAscension.FromDegrees(385 - 360).Radians).Within(1E-9));
    }
    
    [Test]
    public void Property_Tests() {
        var ra1 = new RightAscension(0.875);
        Assert.That(ra1.Radians, Is.EqualTo(0.875));
        Assert.That(ra1.Degrees, Is.EqualTo(0.875 * Constants.DegreesPerRadian).Within(1E-9));
        Assert.That(ra1.Arcminutes, Is.EqualTo(60.0 * 0.875 * Constants.DegreesPerRadian).Within(1E-9));
        Assert.That(ra1.Arcseconds, Is.EqualTo(3600.0 * 0.875 * Constants.DegreesPerRadian).Within(1E-9));
        Assert.That(ra1.Hours, Is.EqualTo(ra1.Degrees / 15.0).Within(1E-9));
        Assert.That(ra1.WrapMin, Is.EqualTo(0));
        Assert.That(ra1.WrapMax, Is.EqualTo(2 * Math.PI));
    }
    
    [Test]
    public void Method_Tests() {
        // .Equals():
        var ra1 = new RightAscension(0.875);
        var ra2 = new RightAscension(0.875);
        var ra3 = new RightAscension(0.111);
        Assert.That(ra1, Is.EqualTo(ra2));
        Assert.That(ra1, Is.Not.EqualTo(ra3));

        // .ToString():
        var words = ra1.ToString().Split();
        Assert.That(words[0], Is.EqualTo("RightAscension:"));
        Assert.That(words[1], Is.EqualTo("0.875"));
        Assert.That(words[2], Is.EqualTo("radians"));
        Assert.That(words[3], Is.EqualTo("="));
        Assert.That(words[4].StartsWith("50.1338"), Is.True);
        Assert.That(words[5].StartsWith("degrees"));
    }
}
