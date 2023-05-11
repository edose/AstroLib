using System.Text.RegularExpressions;

/// @namespace AstroLib.Catalog.Atlas2
/// All classes that directly handle the ATLAS refcat2 star catalog.
namespace AstroLib.Catalog.Atlas2;

/// <summary>Provides user access to an installation of the
/// Atlas refcat2 catalog on local storage. </summary>
public class Atlas2Catalog {

    public string TopPath { get; init; }    
    
    public static readonly DateTime CatalogEpoch = 
        new DateTime(2015, 7, 1, 0, 0, 0);

    public static readonly int RequiredFieldsPerDataLine = 44;

    public static readonly List<string> RecognizedSubdirectoryNames = new() {
        "mag-0-16", "mag-16-17", "mag-17-18",
        "mag-18-19", "mag-19-20" };
    
    public static readonly Dictionary<string, string> SubdirectoryCodes = new() {
        {"mag-0-16", "A"}, {"mag-16-17", "B"}, {"mag-17-18", "C"},
        {"mag-18-19", "D"}, {"mag-19-20", "E"} };

    public static readonly Dictionary<string, int> ExpectedFileCounts = new() {
        {"mag-0-16", 64800}, {"mag-16-17", 64799}, {"mag-17-18", 64800},
        {"mag-18-19", 64800}, {"mag-19-20", 64800} };

    public List<string> CatalogSubdirectoryNamesPresent { get; init; }
    public List<string> CatalogSubdirectoryPathsPresent { get; init; }

    /// <summary>Private default constructor, to prevent user construction without input data.</summary>
    private Atlas2Catalog() {
    }

    /// <summary>Normal constructor.</summary>
    /// <param name="topPath">Path to top of Atlas refcat2 catalog.</param>
    public Atlas2Catalog(string topPath) {
        TopPath = topPath;
        var allPathsPresent = Directory.GetDirectories(TopPath,
            "mag-*", SearchOption.TopDirectoryOnly);
        var allNamesPresent = allPathsPresent
            .Select(s=>s.Split(@"\").Last()).ToList();
        CatalogSubdirectoryPathsPresent = new();
        CatalogSubdirectoryNamesPresent = new();
        // We go to these lengths so that both Catalog...Present lists
        // are in same order as RecognizedSubdirectories.
        foreach (var recognized in RecognizedSubdirectoryNames) {
            var indexFound = allNamesPresent.IndexOf(recognized);
            if (indexFound != -1) {
                CatalogSubdirectoryPathsPresent.Add(allPathsPresent[indexFound]);
                CatalogSubdirectoryNamesPresent.Add(allNamesPresent[indexFound]);
            }
        }
    }

    /// <summary> Runs several checks to raise confidence that this catalog installation
    /// is ready for use.</summary>
    /// <returns>True iff catalog installation appears ready.</returns>
    public List<string> ValidateCatalogInstallation() {
        var warnings = new List<string>();

        // Report absence of catalog subdirectories; in which case, return immediately:
        if (CatalogSubdirectoryNamesPresent.Count < 1) {
            warnings.Add("No catalog subdirectories found.");
            return warnings;
        }

        // Report any subdirectories with invalid number of .rc2 files:
        foreach (var sd in CatalogSubdirectoryNamesPresent) {
            var sdPath = Path.Combine(TopPath, sd);
            var actualFileCount = Directory.EnumerateFiles(sdPath, "*.rc2").Count();
            var expectedFileCount = ExpectedFileCounts[sd];
            if (actualFileCount != expectedFileCount)
                warnings.Add($"In subdirectory {sdPath}: " +
                             $"expected {expectedFileCount} .rc2 files but found {actualFileCount}.");
        }

        // Report any files ending in .rc2 but having invalid filenames (this could be slow):
        const string validFilenamePattern = @"[0123]\d{2}[\+\-]\d{2}.rc";
        foreach (var path in CatalogSubdirectoryPathsPresent) {
            foreach (var fullpath in Directory.EnumerateFiles(path, "*.rc2")) {
                if (Regex.IsMatch(Path.GetFileName(fullpath), validFilenamePattern, RegexOptions.Compiled)) 
                    continue;
                if (fullpath.EndsWith(".rc2")) {
                    warnings.Add($"Invalid .rc2 filename in fullpath {fullpath}.");
                }
            }
        }

        // Verify (cursorily) syntax of 3 .rc2 files in each existing catalog subdirectory:
        foreach (var path in CatalogSubdirectoryPathsPresent) {
            var filenames = Directory.EnumerateFiles(path, "*.rc2").ToArray();
            var filenamesToTest = new string[] {
                filenames.First(), filenames.Last(),
                filenames[filenames.Length / 3]
            };
            foreach (var fn in filenamesToTest) {
                var lines = readLinesFromOneAtlasFile(path, fn);
                foreach (var line in lines) {
                    var fields = line.Split(',');
                    if (fields.Count() != RequiredFieldsPerDataLine) {
                        warnings.Add($"Invalid number of fields per line ({fields.Count()}) " +
                                     $"in subdirectory {path}, filename {fn}.");
                        foreach (var field in fields) {
                            if (!int.TryParse(field, out _)) {
                                warnings.Add(
                                    $"Non-integer field {field} in subdirectory {path}, filename {fn}.");
                            }
                        }
                    }
                }

                
            }
        }
        return warnings;
    }
    
    /// <summary>Returns list of strings from one degree-square file specified by its filename.
    /// Each string in the list is one line read from that file.</summary>
    /// <param name="subdirectory">Atlas refcat2 subdirectory, e.g. "mag-0-16".</param>
    /// <param name="filename">Filename, e.g., "233+04.rc2" </param>
    /// <returns></returns>
    private List<string> readLinesFromOneAtlasFile(string subdirectory, string filename) {
        var fullPath = Path.Combine(TopPath, subdirectory, filename);
        return (File.Exists(fullPath)) ? File.ReadLines(fullPath).ToList() : new List<string>();
    }

    /// <summary>Returns list of strings from one degree-square file specified by RA and Dec (integers).
    /// Each string in the list is one line read from that file.</summary>
    /// <param name="subdirectory">Atlas refcat2 subdirectory.</param>
    /// <param name="ra">Right Ascension integer.</param>
    /// <param name="dec">Declination integer.</param>
    /// <returns></returns>
    public List<string> readLinesFromOneAtlasFile(string subdirectory, int ra, int dec) {
        var filename = $"{ra:000}{dec:+00;-00}.rc2";
        return readLinesFromOneAtlasFile(subdirectory, filename);
    }


    public override string ToString() {
        return $"Atlas2Catalog object from {TopPath}.";
    }
}
