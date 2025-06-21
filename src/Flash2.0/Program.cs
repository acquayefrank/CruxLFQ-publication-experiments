using FlashLFQ;

class Program
{
    static void Main(string[] args)
    {
        if (args.Length < 2)
        {
            Console.WriteLine("Please provide the list of files and the file path as command line arguments.");
            return;
        }

        Console.WriteLine("Processing FlashLFQ!");

        var files = args.Take(args.Length - 1).ToList();
        string filePath = args.Last();

        int fileCounter = 0;
        List<SpectraFileInfo> specFiles = new List<SpectraFileInfo>();
        foreach (var file in files)
        {
            SpectraFileInfo specFile = new SpectraFileInfo(file, "a", fileCounter, 0, 0);
            specFiles.Add(specFile);
            fileCounter++;
        }

        List<Identification> identifications = new List<Identification>();

        try
        {
            using StreamReader reader = new StreamReader(filePath);
            string line;
            int lineNumber = 0;
            while ((line = reader.ReadLine()) != null)
            {
                if (lineNumber != 0)
                {
                    string[] values = line.Split('\t');
                    int counter = 0;
                    string sequence = "", modifiedSequence = "", protein_id = "";
                    double monoisotopicMass = 0;
                    double ms2RetentionTimeInMinutes = 0;
                    int chargeState = 0;
                    bool shouldNotSkip = true;
                    foreach (string value in values)
                    {
                        if (counter == 11)
                        {
                            double q_value = double.Parse(value);
                            if (q_value > 0.01)
                            {
                                shouldNotSkip = false;
                            }
                            else
                            {
                                shouldNotSkip = true;
                            }
                        }
                        if (counter == 17) { sequence = value; modifiedSequence = value; }
                        if (counter == 2) { chargeState = int.Parse(value); }
                        if (counter == 6) { monoisotopicMass = double.Parse(value); }
                        if (counter == 3) { ms2RetentionTimeInMinutes = double.Parse(value); }
                        if (counter == 19)
                        {
                            protein_id = value;
                        }

                        counter++;
                    }
                    if (shouldNotSkip)
                    {
                        var pg = new ProteinGroup("Protein", "gene", protein_id);

                        foreach (SpectraFileInfo specFile in specFiles)
                        {
                            Identification id = new Identification(specFile, sequence, modifiedSequence, monoisotopicMass, ms2RetentionTimeInMinutes, chargeState, new List<ProteinGroup> { pg });
                            identifications.Add(id);
                        }
                    }
                }
                lineNumber++;
            }
        }
        catch (Exception ex)
        {
            Console.WriteLine($"An error occurred while reading the file: {ex.Message}");
        }

        FlashLfqEngine engine = new FlashLfqEngine(identifications, normalize: false, maxThreads: 1);

        var results = engine.Run();
        string numFiles = files.Count.ToString();
        results.WriteResults(
            peaksOutputPath: $"{numFiles}_FlashLFQ+mods+protein_id_peaks.txt",
            modPeptideOutputPath: $"{numFiles}_FlashLFQ+mods+protein_id_modpep.txt",
            $"{numFiles}_FlashLFQ+mods+protein_id_prot.txt",
            $"{numFiles}_FlashLFQ+mods+protein_id_bays.txt",
            false
        );
    }
}