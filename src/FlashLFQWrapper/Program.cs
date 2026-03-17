using FlashLFQ;
using System.IO;
using System.Collections.Generic;
using System.Globalization;
using System.Diagnostics;

class Program
{
    static void Main(string[] args)
    {
        if (args.Length < 1)
        {
            Console.WriteLine("Usage: FlashLFQWrapper <input_file> [--out <output_folder>] [--replicates <replicate_file>]");
            return;
        }

        Console.WriteLine("Processing FlashLFQ!");

        // Parse --out argument at any position
        string outputFolder = "flash_lfq_results";
        List<string> argList = args.ToList();
        int outIndex = argList.IndexOf("--out");
        if (outIndex != -1)
        {
            if (outIndex + 1 >= argList.Count)
            {
                Console.WriteLine("Error: --out flag provided but no folder specified.");
                return;
            }
            outputFolder = argList[outIndex + 1];
            argList.RemoveAt(outIndex + 1);
            argList.RemoveAt(outIndex);
        }

        // Parse --replicates argument
        string replicatesFile = null;
        int repIndex = argList.IndexOf("--replicates");
        if (repIndex != -1)
        {
            if (repIndex + 1 >= argList.Count)
            {
                Console.WriteLine("Error: --replicates flag provided but no file specified.");
                return;
            }
            replicatesFile = argList[repIndex + 1];
            argList.RemoveAt(repIndex + 1);
            argList.RemoveAt(repIndex);
        }

        if (argList.Count < 1)
        {
            Console.WriteLine("Error: Not enough arguments after parsing --out.");
            return;
        }

        string filePath = argList.Last();

        // Ensure output folder exists
        if (!Directory.Exists(outputFolder))
        {
            Directory.CreateDirectory(outputFolder);
        }

        // Load replicate metadata from file if provided
        var replicateMap = new Dictionary<string, (string condition, int biologicalReplicate, int fraction, int technicalReplicate)>();
        if (replicatesFile != null)
        {
            foreach (var repLine in File.ReadLines(replicatesFile).Skip(1))
            {
                var cols = repLine.Split('\t');
                if (cols.Length < 5) continue;
                string path = cols[0].Trim();
                string cond = cols[1].Trim();
                if (!int.TryParse(cols[2].Trim(), out int bioRep)) continue;
                if (!int.TryParse(cols[3].Trim(), out int frac)) continue;
                if (!int.TryParse(cols[4].Trim(), out int techRep)) continue;
                replicateMap[path] = (cond, bioRep, frac, techRep);
            }
        }

        // Read all lines once (except header)
        var allLines = File.ReadLines(filePath).Skip(1).ToList();

        // Collect unique spectral files
        var spectralFiles = new HashSet<string>(
            allLines.Select(line => line.Split('\t'))
                    .Where(values => values.Length >= 10)
                    .Select(values => values[5])
        );

        // Create SpectraFileInfo objects from replicate file or regex fallback
        var specFileDict = spectralFiles.ToDictionary(
            sf => sf,
            sf =>
            {
                if (replicateMap.TryGetValue(sf, out var meta))
                {
                    return new SpectraFileInfo(sf, meta.condition, meta.biologicalReplicate, meta.fraction, meta.technicalReplicate);
                }

                // Fallback: extract condition and replicate from filename
                string fileName = Path.GetFileName(sf);
                string condition = "a";
                int biologicalReplicate = 0;

                var matches = System.Text.RegularExpressions.Regex.Matches(fileName, @"([A-E])(\d+)");
                if (matches.Count > 0)
                {
                    var match = matches[matches.Count - 1];
                    condition = match.Groups[1].Value;
                    if (int.TryParse(match.Groups[2].Value, out int replicate))
                        biologicalReplicate = replicate - 1;
                }

                return new SpectraFileInfo(sf, condition, biologicalReplicate, 0, 0);
            }
        );

        // Create identifications
        var identifications = new List<Identification>();
        foreach (var line in allLines)
        {
            var values = line.Split('\t');
            if (values.Length < 10) continue;

            if (!double.TryParse(values[2], NumberStyles.Any, CultureInfo.InvariantCulture, out double monoIsotopicMass)) continue;
            if (!int.TryParse(values[4], out int precursorCharge)) continue;
            if (!double.TryParse(values[6], NumberStyles.Any, CultureInfo.InvariantCulture, out double ms2RetentionTimeInMinutes)) continue;

            string sequence = values[0];
            string spectralFile = values[5];
            string modifications = values[8];
            string protein_id = values[9];

            if (!specFileDict.TryGetValue(spectralFile, out var specFile)) continue;

            var pg = new ProteinGroup("Protein", "gene", protein_id);
            var id = new Identification(
                specFile,
                sequence,
                modifications,
                monoIsotopicMass,
                ms2RetentionTimeInMinutes,
                precursorCharge,
                new List<ProteinGroup> { pg }
            );
            identifications.Add(id);
        }

        long memoryBefore = Process.GetCurrentProcess().PrivateMemorySize64;

        FlashLfqEngine engine = new FlashLfqEngine(identifications, normalize: true, maxThreads: 1);

        long memoryAfter = Process.GetCurrentProcess().PrivateMemorySize64;
        long memoryUsed = memoryAfter - memoryBefore;
        Console.WriteLine($"Approx. memory used: {memoryUsed / 1024.0:F2} KB");

        var results = engine.Run();
        results.WriteResults(
            peaksOutputPath: Path.Combine(outputFolder, "FlashLFQ+mods+protein_id_peaks.txt"),
            modPeptideOutputPath: Path.Combine(outputFolder, "FlashLFQ+mods+protein_id_modpep.txt"),
            Path.Combine(outputFolder, "FlashLFQ+mods+protein_id_prot.txt"),
            Path.Combine(outputFolder, "FlashLFQ+mods+protein_id_bays.txt"),
            false
        );
    }
}