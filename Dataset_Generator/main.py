from Bio import Entrez, SeqIO
from Bio.Blast import NCBIWWW, NCBIXML

class DatasetGenerator:
    """
    This class generates fasta files of n sequences. The sequences are grouped into n_cl clusters.
    Each cluster is defined by a seed sequence and contains very similar sequences.
    The remaining sequences of a cluster are obtained by blasting the seed sequence.
    """    
    
    def generate(self, n, n_cl, seeds):
        """
        Generates a fasta file with n sequences grouped into n_cl clusters.
        :param n: number of sequences
        :param n_cl: number of clusters
        :param seeds: list of seed sequences (uids)
        :return: path to the generated fasta file
        """
        self._check_for_correct_input(n, n_cl, seeds)
        fasta_content = ""
        for i in range(n_cl):
            records_of_cluster_i = []
            seed_uid = seeds[i]
            # download seed sequence
            seed_record = self._download_sequences([seed_uid])
            records_of_cluster_i.append(seed_record)
            # blast for n/n_cl - 1 similar sequences
            print(f"Finding similar sequences for cluster {i + 1}...")
            records_of_similar_sequences = self.get_records_of_similar_sequences(seed_record, n // n_cl - 1)
            records_of_cluster_i.extend(records_of_similar_sequences)
            # add sequences of cluster i to fasta_content and mark the cluster
            for record in records_of_cluster_i:
                fasta_content += f">Cluster_{i+1}_Sequence_{record.id}\n{str(record.seq)}\n"
        # save the generated fasta content to a file
        output_path = f"example_files/generated_fasta/{n}_seq_{n_cl}_clusters.fasta"
        print(f"Saving fasta file at {output_path}")
        self._save_to_file(output_path, fasta_content)
        return output_path
    
    
    def _check_for_correct_input(self, n, n_cl, seeds):
        """
        Checks if the input parameters are valid.
        :param n: number of sequences
        :param n_cl: number of clusters
        :param seeds: list of seed sequences (uids)
        """
        if n <= 0:
            raise ValueError("Number of sequences must be positive.")
        if n_cl > n:
            raise ValueError("Number of clusters cannot be greater than number of sequences.")
        if len(seeds) != n_cl:
            raise ValueError("Number of seed sequences must match number of clusters.")
        
    
    def _download_sequences(self, uids: list):
        """
        Downloads a protein sequence from NCBI given its UID.
        :param uid: the UID of the sequence
        :return: the sequence as a string
        """
        Entrez.email = "aron.wichtner@tuebingen.mpg.de"
        uids_as_string = ",".join(uids)
        handle = Entrez.efetch(db="protein", id=uids_as_string, rettype="fasta", retmode="text")
        records = list(SeqIO.parse(handle, "fasta"))
        handle.close()
        return records
    
    
    def _blast_sequence(self, record, size_of_hitlist):
        """
        PSI-Blasts a protein record.
        :param record: the protein record to blast
        :return: result of BLAST
        """
        # Run PSI-BLAST with the fetched sequence
        result_handle = NCBIWWW.qblast(
            program="blastp",
            database="nr",
            sequence=str(record.seq),
            service="psi",
            hitlist_size=size_of_hitlist,
            word_size=3
        )
        return result_handle
    
    
    def get_records_of_similar_sequences(self, record, num_sequences):
        """
        Retrieves records of sequences similar to the given record sequence.
        :param record: the record to compare
        :param num_sequences: number of similar records to retrieve
        :return: list of SeqRecord objects
        """
        results_handle = self._blast_sequence(record, num_sequences)
        self.psiblast_results_path = "Dataset_Generator/psiblast_results.xml"
        self._save_to_file(self.psiblast_results_path, results_handle.read())
        results_handle.close()
        similar_records = self._parse_hits()
        return similar_records
        
        
    def _parse_hits(self):
        """
        Parses the BLAST hits from the XML result file.
        :return: list of SeqRecord objects
        """
        fasta_records = []
        hit_ids = []
        with open(self.psiblast_results_path) as f:
            blast_records = NCBIXML.parse(f)
            for blast_record in blast_records:
                # Collect hit IDs
                for alignment in blast_record.alignments:
                    hit_ids.append(alignment.hit_id)
            fasta_records = self._download_sequences(hit_ids)
        return fasta_records


    def _save_to_file(self, filename, content):
        """
        Saves a file with a given content.
        :param filename: path to the file
        :param content: content to write into the file
        """
        with open(filename, "w") as f:
            f.write(content)
            
            
# test
#"Q99895", "P42212"
seeds = ["P68871"]
generator = DatasetGenerator()
generator.generate(9, 1, seeds)
        
    