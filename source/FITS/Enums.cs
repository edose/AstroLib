namespace AstroLib.Fits;

public enum FitsDataType
// Values equal those encoded in header's 'BITPIX' field, as defined in FITS standard.
{
    Character = 8,
    SignedInt16 = 16,
    SignedInt32 = 32,
    SignedInt64 = 64,
    SingleFloat = -32,
    DoubleFloat = -64
}

public enum ValueType
// Specifies which value type a header record contains.
{
    CharacterString,
    Logical,
    Integer,
    Floating,
    Unknown
}