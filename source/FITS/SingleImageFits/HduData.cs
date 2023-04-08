// ReSharper disable StringLiteralTypo

using System.ComponentModel;
// ReSharper disable IdentifierTypo

namespace AstroLib.Fits.SingleImageFits;

/// <summary>Represents image data of this Single Image FITS file.</summary>
public class HduData
{
    public float[]  DataAsSingleFloat { get; private set; }
    public double[] DataAsDoubleFloat { get; private set; }

    /// <summary>Constructor. Reads, parses, stores image data for this SIF file.</summary>
    public HduData(FileStream fs, HduHeader hduHeader) {
        var naxis1 = int.Abs((int) hduHeader.ValueRead["NAXIS1"].ValueInteger!);
        var naxis2 = int.Abs((int) hduHeader.ValueRead["NAXIS2"].ValueInteger!);
        var pixelCount = naxis1 * naxis2;
        var bitpix = (int) hduHeader.ValueRead["BITPIX"].ValueInteger!;
        var fitsDataType = (FitsDataType) bitpix;
        var bytesPerPixel = int.Abs(bitpix) / 8;
        var byteCount = pixelCount * bytesPerPixel;
        var rawArray = LoadRawArray(fs, byteCount);
        fs.Close();
        
        byte[] bytesLe;
        if (BitConverter.IsLittleEndian)  
            bytesLe = ReverseBytes(rawArray, fitsDataType, bytesPerPixel);
        else
            bytesLe = rawArray;

        var bscale = (double) hduHeader.ValueRead["BSCALE"].ValueFloating!;
        var bzero = (double) hduHeader.ValueRead["BZERO"].ValueFloating!;
        DataAsDoubleFloat = CastDoubleFloatAndScale(bytesLe, fitsDataType, bscale, bzero);
        DataAsSingleFloat = RecastSingleFloat(DataAsDoubleFloat);
    }

    private static byte[] LoadRawArray(FileStream fs, int byteCount) {
        var rawArray = new byte[byteCount];
        var bytesRead = fs.Read(rawArray, 0, byteCount);
        if (bytesRead != byteCount)
            throw new InvalidDataException($"{fs.Name}: {byteCount} data bytes expected, " +
                                           $"but only {bytesRead} could be read.");
        return rawArray;
    }

    private static byte[] ReverseBytes(byte[] bytesAsRead, FitsDataType fitsDataType, int bytesPerValue) {
        // TODO: extract this as a more general fn (not needing FitsDataType).
        // TODO: check for integer number of values in bytesAsRead (or maybe just get valueCount).
        // Values from FITS file are expressed as big-endian per FITS standard.
        // This reverses byte order within values to little-endian,
        //    for x64 and other little-endian environments.
        var pOptions = new System.Threading.Tasks.ParallelOptions {
            MaxDegreeOfParallelism = Math.Max(Environment.ProcessorCount - 1, 1)
        };
        var bytesReversed = new byte[bytesAsRead.Length];
        var valueCount = bytesAsRead.Length / bytesPerValue;        
        switch (fitsDataType) {
            case FitsDataType.Character:
                // one byte per value: no swapping warranted. Rare.
                bytesAsRead.CopyTo(bytesReversed, 0);
                break;
            case FitsDataType.SignedInt16:
                // two bytes per value.
                Parallel.For(0, valueCount, pOptions, iValue => {
                    var i = bytesPerValue * iValue;
                    bytesReversed[i    ] = bytesAsRead[i + 1];
                    bytesReversed[i + 1] = bytesAsRead[i    ];
                });
                break;
            case FitsDataType.SignedInt32:
            case FitsDataType.SingleFloat:
                // four bytes per value.
                Parallel.For(0, valueCount, pOptions, iValue => {
                    var i = bytesPerValue * iValue;
                    bytesReversed[i]     = bytesAsRead[i + 3];
                    bytesReversed[i + 1] = bytesAsRead[i + 2];
                    bytesReversed[i + 2] = bytesAsRead[i + 1];
                    bytesReversed[i + 3] = bytesAsRead[i];
                });
                break;
            case FitsDataType.SignedInt64: // rare.
            case FitsDataType.DoubleFloat:
                // eight bytes per value.
                Parallel.For(0, valueCount, pOptions, iValue => {
                    var i = bytesPerValue * iValue;
                    bytesReversed[i]     = bytesAsRead[i + 7];
                    bytesReversed[i + 1] = bytesAsRead[i + 6];
                    bytesReversed[i + 2] = bytesAsRead[i + 5];
                    bytesReversed[i + 3] = bytesAsRead[i + 4];
                    bytesReversed[i + 4] = bytesAsRead[i + 3];
                    bytesReversed[i + 5] = bytesAsRead[i + 2];
                    bytesReversed[i + 6] = bytesAsRead[i + 1];
                    bytesReversed[i + 7] = bytesAsRead[i];
                });
                break;
            default:
                throw new InvalidEnumArgumentException("Fits.HduData.ReverseBytes()");
        }
        return bytesReversed;
    }

    private static double[] CastDoubleFloatAndScale(byte[] bytesLe, FitsDataType fitsDataType, 
        double bscale, double bzero) {
        return new double[1];  // TODO: replace with actual code.
    }

    private static float[] RecastSingleFloat(double[] doubleArray) {
        var floatArray = new float[doubleArray.Length];
        return new float[1];  // TODO: replace with actual code.
    }

}