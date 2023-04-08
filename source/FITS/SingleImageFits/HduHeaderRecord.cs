using System.Runtime.InteropServices.ComTypes;

namespace AstroLib.Fits.SingleImageFits;

/// <summary>Represents one valid HDU header record.</summary>
public class HduHeaderRecord {
    public string RecordString { get; }
    public string Keyword { get; }
    public bool HasValueIndicator { get; }
    public ValueType ValueType { get; private set; }
    public string? ValueString { get; private set; }
    public bool? ValueLogical { get; private set; }
    public int? ValueInteger { get; private set; }
    public double? ValueFloating { get; private set; }
    public string? Comment { get; set; }

    /// <summary>Primary constructor.</summary>
    /// <param name="recordString">Raw record string of length 80, lifted directly from HDU header.</param>
    public HduHeaderRecord(string recordString) {
        RecordString = recordString;
        Keyword = RecordString.Substring(0, 8).TrimEnd();
        var valueIndicator = RecordString.Substring(8, 2);
        HasValueIndicator = (valueIndicator == "= ");
        ValueType = ValueType.Unknown;
        ValueString = null;
        ValueLogical = null;
        ValueInteger = null;        
        ValueFloating = null;
        Comment = null;

        switch (Keyword) {
            case "CONTINUE":
                (ValueString, Comment) = ParseContinueRecord(RecordString);
                break;
            case "COMMENT":
            case "HISTORY":
            case "":
                (ValueString, Comment) = ParseCommentaryRecord(RecordString);
                break;
            default:
                (var valueType, Comment) = ParseValueKeywordRecord(RecordString);
                break;
        }
    }
    
    /// <summary> Parse a header record with CONTINUE keyword.</summary>
    private static (string?, string?) ParseContinueRecord(string recordString) {
        string? valueString, commentString;
        var valueType = ValueType.Unknown;
        int iRemaining;
       (valueString, valueType, iRemaining) = TryParseCharacterStringValue(recordString, 9);
        commentString = (valueType == ValueType.CharacterString)
            ? ParseComment(recordString, iRemaining)
            : null;
        return (valueString, commentString);
    }

    /// <summary> Parse a "commentary keyword "FITS header record.</summary>
    private static (string?, string?) ParseCommentaryRecord(string recordString) {
        string? valueString = null; // by definition, for commentary record, per FITS standard.
        string? commentString = recordString.Substring(8).TrimEnd();
        if (commentString.Length == 0)
            commentString = null;
        return (valueString, commentString);
    }
    
    /// <summary> Parse a "single record keyword" FITS header record.
    /// Inserts value into class's correct ValueString, ValueLogical, ... location directly (not as return). 
    /// </summary>
    /// <returns>2-Tuple: (valueType, commentString)</returns>
    private (ValueType, string?) ParseValueKeywordRecord(string recordString) {
        string? commentString;
        ValueType valueType;

        (string? stringValue, valueType, var iRemaining) = 
            TryParseCharacterStringValue(recordString, 10);
        if (valueType == ValueType.CharacterString) {
            ValueType = valueType;
            ValueString = stringValue;
            commentString = ParseComment(recordString, iRemaining);
            return (valueType, commentString);
        }
        (bool? logicalValue, valueType, iRemaining) = TryParseLogicalValue(recordString);
        if (valueType == ValueType.Logical) {
            ValueType = valueType;
            ValueLogical = logicalValue;
            commentString = ParseComment(recordString, iRemaining);
            return (valueType, commentString);
        }
        (int? integerValue, valueType, iRemaining) = TryParseIntegerValue(recordString);
        if (valueType == ValueType.Integer) {
            ValueType = valueType;
            ValueInteger = integerValue;
            commentString = ParseComment(recordString, iRemaining);
            return (valueType, commentString);
        }
        (double? floatingValue, valueType, iRemaining) = TryParseFloatingValue(recordString);
        if (valueType == ValueType.Floating) {
            ValueType = valueType;
            ValueFloating = floatingValue;
            commentString = ParseComment(recordString, iRemaining);
            return (valueType, commentString);
        }
        
        // TODO: Consider throwing custom exception here (cannot parse this header record).
        return (ValueType.Unknown, null); // should never get here. 
    }
    
    
    /// <summary>Returns string value for this FITS header record, if value is a character string.</summary>
    /// <returns>3-Tuple: (string? intValue, ValueType valueType, int iRemaining)</returns>
    private static (string?, ValueType, int) TryParseCharacterStringValue(string recordString, int iStart) {
        (string?, ValueType, int) notCharacterStringValue = (null, ValueType.Unknown, iStart);
        
        // Try to find opening single-quote:
        var iOpenQuote = recordString.IndexOf('\'', iStart);
        if (iOpenQuote == -1)
            return notCharacterStringValue;

        // Try to find closing single-quote that is not escaped as a literal quote:
        var iStartCloseQuoteSearch = iOpenQuote + 1;
        while (iStartCloseQuoteSearch < recordString.Length - 1) {
            var iCloseQuote = recordString.IndexOf('\'', iStartCloseQuoteSearch);
            if (iCloseQuote == -1)
                return notCharacterStringValue;
            var previousChar = recordString[iCloseQuote - 1];
            var closeQuoteFound =  (previousChar != '\'') || (iCloseQuote == iOpenQuote + 1);
            if (!closeQuoteFound)
                continue; 
            
            // Extract the character string.
            var rawCharacterStringValue = recordString[(iOpenQuote + 1)..iCloseQuote].TrimEnd(); 
            
            // Collapse all instances of 2 successive quotes to a single quote, per FITS standard.
            var characterStringValue = rawCharacterStringValue.Replace("\'\'", "\'");
            return (characterStringValue, ValueType.CharacterString, iCloseQuote + 1); // normal case.
        }
        return notCharacterStringValue;  // This is not a Character String-valued FITS header record.
    }

    /// <summary>Returns logical/boolean value for this FITS header record, if value is logical.</summary>
    /// <returns>3-Tuple: (bool? logicalValue, ValueType valueType, int iRemaining)</returns>
    private static (bool?, ValueType, int) TryParseLogicalValue(string recordString) {
        var (firstWord, iRemaining) = FindFirstWordInString(recordString, 10);
        if (firstWord == null)
            return (null, ValueType.Unknown, 10);
        var valueString = firstWord.Split('/')[0];
        var isLogical = firstWord is "T" or "F";
        var logicalValue = (firstWord == "T");
        if (isLogical)
            return (logicalValue, ValueType.Logical, iRemaining); // Success.
        return (null, ValueType.Unknown, 10);  // This is not a Logical-value FITS header record.
    }

    /// <summary>Returns integer value for this FITS header record, if value is integer.</summary>
    /// <returns>3-Tuple: (int? intValue, ValueType valueType, int iRemaining)</returns>
    private static (int?, ValueType, int) TryParseIntegerValue(string recordString) {
        var (firstWord, iRemaining) = FindFirstWordInString(recordString, 10);
        if (firstWord == null)
            return (null, ValueType.Unknown, 10);
        var valueString = firstWord.Split('/')[0];
        var isInteger = int.TryParse(valueString, out var intValue);
        if (isInteger)
            return (intValue, ValueType.Integer, iRemaining); // Success.
        return (null, ValueType.Unknown, 10);  // This is not an Integer-value FITS header record.
    }

    /// <summary>Returns floating(-point) value for this FITS header record,
    /// if value is floating(-point).</summary>
    /// <returns>3-Tuple: (double? intValue, ValueType valueType, int iRemaining)</returns>
    private static (double?, ValueType, int) TryParseFloatingValue(string recordString) {
        var (firstWord, iRemaining) = FindFirstWordInString(recordString, 10);
        if (firstWord == null)
            return (null, ValueType.Unknown, 10);
        var valueString = firstWord.Split('/')[0];
        var isFloating = double.TryParse(valueString, out var floatingValue);
        if (isFloating)
            return (floatingValue, ValueType.Floating, iRemaining); // Success.
        return (null, ValueType.Unknown, 10);  // This is not a Floating-value FITS header record.
    }
    
    private static (string?, int) FindFirstWordInString(string recordSubstring, int iStart) {
        if (iStart > recordSubstring.Length - 1 || iStart < 0)
            return (null, iStart);
        var substringToSearch = recordSubstring.Substring(iStart);
        var rawFirstWord = substringToSearch.TrimStart().Split(' ')[0];
        var trimmedFirstWord = rawFirstWord.Trim();
        if (trimmedFirstWord == "")
            trimmedFirstWord = null;
        trimmedFirstWord = (trimmedFirstWord != "") ? trimmedFirstWord : null;
        var iRemaining = iStart + 
                            substringToSearch.IndexOf(rawFirstWord, StringComparison.Ordinal) + 
                            rawFirstWord.Length;
        return (trimmedFirstWord, iRemaining);
    }
    
    private static string? ParseComment(string recordSubstring, int iStart) {
        if (iStart > recordSubstring.Length - 1)
            return null;
        var possibleSubstring = recordSubstring.Substring(iStart);
        var iSlash = possibleSubstring.IndexOf('/');
        var commentAvailable = ((iSlash != -1) && (iSlash+1 < possibleSubstring.Length - 1));
        var commentString = commentAvailable ? possibleSubstring.Substring(iSlash + 1).TrimEnd() : null;
        return commentString;
    }


    
    
}

