namespace AstroLib.Fits.SingleImageFits;

/// <summary>Container for header of (sole) HDU of this Single Image FITS file.</summary>
public class HduHeader
{
    private const int BytesPerHeaderBlock = 2880;
    private const int BytesPerHeaderRecord = 80;
    private const int RecordsPerHeaderBlock = BytesPerHeaderBlock / BytesPerHeaderRecord; // 36

    public List<HduHeaderRecord> RecordList { get; private set; }
    public Dictionary<string, HduHeaderRecord> ValueRead { get; private set; }
    public bool MandatoryKeywordsPresent { get; private set; }
    // public Dictionary<string, string?> Values { get; private set; }

    /// <summary>Constructor. Reads, parses, stores header records for this Single Image FITS file.
    /// Verifies that header info specifies a Single Image FITS.</summary>
    public HduHeader(FileStream fs)
    {
        // Temporary backing fields for read-only Properties.
        int bytesConsumed = 0;  // includes entire header blocks, not just portions resulting in records.
        var hduHeaderRecords = new List<HduHeaderRecord>();
        var hduHeaderDictionary = new Dictionary<string, HduHeaderRecord>();
        var hduHeaderValues = new Dictionary<string, string?>();

        var hduHeaderRecordStrings = ExtractHeaderStrings(fs);

        // Parse raw header strings into a list of header records:
        foreach (var rs in hduHeaderRecordStrings)
            hduHeaderRecords.Add(new HduHeaderRecord(rs));

        // Make dictionaries for (1) HDU records, and (2) HDU record value strings:
        hduHeaderDictionary = MakeHduHeaderDictionary(hduHeaderRecords);
        var mandatoryKeywordsPresent = VerifyMandatoryKeywords(hduHeaderDictionary);
        // foreach (var record in hduHeaderDictionary)
        //     hduHeaderValues.Add(record.Key, record.Value.ValueString);

        // Fill class properties:
        RecordList = hduHeaderRecords;
        ValueRead = hduHeaderDictionary;
        MandatoryKeywordsPresent = mandatoryKeywordsPresent;
    }

    /// <summary>Chop header into raw record strings, return list of them.</summary>
    public static List<string> ExtractHeaderStrings(FileStream fs)
    {
        Boolean endKeywordIsFound = false;
        var rawEndKeyword = "END".PadRight(8);
        var hduHeaderRecordStrings = new List<string>();
        var headerBlock = new byte[BytesPerHeaderBlock];

        while (!endKeywordIsFound)
        {
            var blockBytesRead = fs.Read(headerBlock, offset: 0, count: BytesPerHeaderBlock);
            if (blockBytesRead != BytesPerHeaderBlock)
                throw new FileLoadException($"Defective Header block--its length " +
                                            $"should be {BytesPerHeaderBlock}");
            for (var i = 0; i < RecordsPerHeaderBlock; i++)
            {
                var recordString = System.Text.Encoding.ASCII.GetString(headerBlock,
                    i * BytesPerHeaderRecord, BytesPerHeaderRecord);
                if (recordString.Trim() != "") hduHeaderRecordStrings.Add(recordString);
                if (!recordString.StartsWith(rawEndKeyword)) continue;
                endKeywordIsFound = true;
                break;
            }
        }
        return hduHeaderRecordStrings;
    }

    /// <summary>Return dictionary (key, record) of header records other than COMMENT & CONTINUE.
    /// Adds only first occurrence of each keyword (to preserve Mandatory Keyword records).</summary>
    public static Dictionary<string, HduHeaderRecord> MakeHduHeaderDictionary
        (List<HduHeaderRecord> hduHeaderRecords)
    {
        var hduHeaderDictionary = new Dictionary<string, HduHeaderRecord>();
        foreach (var r in hduHeaderRecords)
            if (r.Keyword != "COMMENT" && r.Keyword != "CONTINUE")
                hduHeaderDictionary.TryAdd(r.Keyword, r);
        return hduHeaderDictionary;
    }

    /// <summary>Verify presence of all Mandatory Records (per FITS standard)</summary>
    public static bool VerifyMandatoryKeywords(Dictionary<string, HduHeaderRecord> hduHeaderDict)
    {
        string[] mandatoryKeywords = {"SIMPLE", "BITPIX", "NAXIS", "NAXIS1", "NAXIS2", "END"};
        var absentMandatoryKeywords = new List<string>();
        foreach (var mr in mandatoryKeywords)
            if (!hduHeaderDict.ContainsKey(mr))
                absentMandatoryKeywords.Add(mr);
        var mandatoryKeywordsPresent = (absentMandatoryKeywords.Count == 0);
        return mandatoryKeywordsPresent;
        // if (absentMandatoryKeywords.Count > 0)
        //     throw new FileLoadException($"FITS mandatory keywords missing: " +
        //                                 $"{string.Join(" ", absentMandatoryKeywords)}");
    }
}