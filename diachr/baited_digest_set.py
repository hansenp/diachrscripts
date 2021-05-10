from .baited_digest import BaitedDigest
from .diachromatic_interaction import DiachromaticInteraction11
from .diachromatic_interaction_set import DiachromaticInteractionSet
from collections import defaultdict

class BaitedDigestSet:

    def __init__(self):

        # Dictionary that contains all BaitedDigest objects
        self._baited_digest_dict = defaultdict(BaitedDigest)

        # Dictionary with infos about ingestion of interactions
        self._ingest_interaction_set_info_dict = {
            'TOTAL_INTERACTIONS_READ': 0,
            'DI': {
                'NN': 0,
                'NE': 0,
                'EN': 0,
                'EE': 0
            },
            'UIR': {
                'NN': 0,
                'NE': 0,
                'EN': 0,
                'EE': 0
            },
            'UI': {
                'NN': 0,
                'NE': 0,
                'EN': 0,
                'EE': 0
            },
            'ALL': {
                'NN': 0,
                'NE': 0,
                'EN': 0,
                'EE': 0
            },
            'BAITED_DIGESTS': 0
        }

    def ingest_interaction_set(self, d11_inter_set: DiachromaticInteractionSet, verbose: bool = False):

        if verbose:
            print("[INFO] Reading interactions and group them according to baited digests ...")

        for d11_inter in d11_inter_set.interaction_list:
            self._ingest_interaction_set_info_dict['TOTAL_INTERACTIONS_READ'] += 1
            if verbose:
                if self._ingest_interaction_set_info_dict['TOTAL_INTERACTIONS_READ'] % 1000000 == 0:
                    print("\t[INFO] Read " + "{:,}".format(self._ingest_interaction_set_info_dict['TOTAL_INTERACTIONS_READ']) + " interactions ...")

            # Count interaction type
            self._ingest_interaction_set_info_dict[d11_inter.get_category()][d11_inter.enrichment_status_tag_pair] += 1
            self._ingest_interaction_set_info_dict['ALL'][d11_inter.enrichment_status_tag_pair] += 1

            # Skip NN and EE interactions
            if d11_inter.enrichment_status_tag_pair == 'NN' or d11_inter.enrichment_status_tag_pair == 'EE':
                continue

            # Get key from interaction
            if d11_inter.enrichment_status_tag_pair == 'EN':
                enriched_digest_coords = d11_inter.chrA + '\t' + str(d11_inter.fromA) + '\t' + str(d11_inter.toA)
            if d11_inter.enrichment_status_tag_pair == 'NE':
                enriched_digest_coords = d11_inter.chrB + '\t' + str(d11_inter.fromB) + '\t' + str(d11_inter.toB)

            # Create new BaitedDigest object if it doesn't already exists
            if enriched_digest_coords not in self._baited_digest_dict.keys():
                self._baited_digest_dict[enriched_digest_coords] = BaitedDigest()
                self._ingest_interaction_set_info_dict['BAITED_DIGESTS'] += 1

            # Add interaction BaitedDigest object
            self._baited_digest_dict[enriched_digest_coords].add_interaction(d11_inter)

        if verbose:
            print("\t[INFO] Total number of interactions read: " + "{:,}".format(self._ingest_interaction_set_info_dict['TOTAL_INTERACTIONS_READ']))
            print("\t[INFO] Total number of baited digests: " + "{:,}".format(self._ingest_interaction_set_info_dict['BAITED_DIGESTS']))
            print("[INFO] ... done.")

    def get_ingest_interaction_set_info_dict(self):
        """
        :return: Dictionary that contains information about parsed interaction data.
        """
        return self._ingest_interaction_set_info_dict

    def get_ingest_interaction_set_info_report(self):

        report = "[INFO] Report on ingestion of interactions:" + '\n'
        report += "\t[INFO] Total number of interactions read: " + "{:,}".format(self._ingest_interaction_set_info_dict['TOTAL_INTERACTIONS_READ']) + '\n'
        discarded = self._ingest_interaction_set_info_dict['TOTAL_INTERACTIONS_READ'] - \
                    self._ingest_interaction_set_info_dict['ALL']['NE'] - \
                    self._ingest_interaction_set_info_dict['ALL']['EN']
        report += "\t[INFO] Discarded NN and EE interactions: " + "{:,}".format(discarded) + '\n'
        ingested = self._ingest_interaction_set_info_dict['ALL']['NE'] + \
                   self._ingest_interaction_set_info_dict['ALL']['EN']
        report += "\t[INFO] Total number of ingested NE and EN interactions: " + "{:,}".format(ingested) + '\n'
        report += "\t[INFO] Broken down by interaction category and enrichment status: " + '\n'
        report += "\t\t[INFO] DI: " + '\n'
        report += "\t\t\t[INFO] NE: " + "{:,}".format(self._ingest_interaction_set_info_dict['DI']['NE']) + '\n'
        report += "\t\t\t[INFO] EN: " + "{:,}".format(self._ingest_interaction_set_info_dict['DI']['EN']) + '\n'
        report += "\t\t[INFO] UIR: " + '\n'
        report += "\t\t\t[INFO] NE: " + "{:,}".format(self._ingest_interaction_set_info_dict['UIR']['NE']) + '\n'
        report += "\t\t\t[INFO] EN: " + "{:,}".format(self._ingest_interaction_set_info_dict['UIR']['EN']) + '\n'
        report += "\t\t[INFO] UI: " + '\n'
        report += "\t\t\t[INFO] NE: " + "{:,}".format(self._ingest_interaction_set_info_dict['UI']['NE']) + '\n'
        report += "\t\t\t[INFO] EN: " + "{:,}".format(self._ingest_interaction_set_info_dict['UI']['EN']) + '\n'
        report += "\t\t[INFO] ALL: " + '\n'
        report += "\t\t\t[INFO] NE: " + "{:,}".format(self._ingest_interaction_set_info_dict['ALL']['NE']) + '\n'
        report += "\t\t\t[INFO] EN: " + "{:,}".format(self._ingest_interaction_set_info_dict['ALL']['EN']) + '\n'
        report += "\t[INFO] Total number of baited digests: " + "{:,}".format(self._ingest_interaction_set_info_dict['BAITED_DIGESTS']) + '\n'
        report += "[INFO] End of report." + '\n'

        return report

    def get_ingest_interaction_set_table_row(self):

        # Header row
        table_row = ":TR_INGESTION:" + '\t'

        table_row += "TOTAL_INTERACTIONS_READ" + '\t'
        table_row += "DISCARDED" + '\t'
        table_row += "INGESTED" + '\t'

        table_row += "DI_NE" + '\t'
        table_row += "DI_EN" + '\t'

        table_row += "UIR_NE" + '\t'
        table_row += "UIR_EN" + '\t'

        table_row += "UI_NE" + '\t'
        table_row += "UI_EN" + '\t'

        table_row += "ALL_NE" + '\t'
        table_row += "ALL_EN" + '\t'

        table_row += "BAITED_DIGESTS" + '\n'

        # Row with values
        table_row += ":TR_INGESTION:" + '\t'

        table_row += str(self._ingest_interaction_set_info_dict['TOTAL_INTERACTIONS_READ']) + '\t'
        discarded = self._ingest_interaction_set_info_dict['TOTAL_INTERACTIONS_READ'] - \
                    self._ingest_interaction_set_info_dict['ALL']['NE'] - \
                    self._ingest_interaction_set_info_dict['ALL']['EN']
        table_row += str(discarded) + '\t'
        ingested = self._ingest_interaction_set_info_dict['ALL']['NE'] + \
                   self._ingest_interaction_set_info_dict['ALL']['EN']
        table_row += str(ingested) + '\t'

        table_row += str(self._ingest_interaction_set_info_dict['DI']['NE']) + '\t'
        table_row += str(self._ingest_interaction_set_info_dict['DI']['EN']) + '\t'

        table_row += str(self._ingest_interaction_set_info_dict['UIR']['NE']) + '\t'
        table_row += str(self._ingest_interaction_set_info_dict['UIR']['EN']) + '\t'

        table_row += str(self._ingest_interaction_set_info_dict['UI']['NE']) + '\t'
        table_row += str(self._ingest_interaction_set_info_dict['UI']['EN']) + '\t'

        table_row += str(self._ingest_interaction_set_info_dict['ALL']['NE']) + '\t'
        table_row += str(self._ingest_interaction_set_info_dict['ALL']['EN']) + '\t'

        table_row += str(self._ingest_interaction_set_info_dict['BAITED_DIGESTS']) + '\n'

        return table_row

    def proportion_of_directed_interactions_on_individual_chromosomes(self):

        baited_digest_sets_chr = {}

        # Divide BaitedDigestSet according to chromosomes
        for key in self._baited_digest_dict.keys():
            chr, sta, end = key.split('\t')
            if chr not in baited_digest_sets_chr.keys():
                baited_digest_sets_chr[chr] = [self._baited_digest_dict[key]]
            else:
                baited_digest_sets_chr[chr].append(self._baited_digest_dict[key])

        for chr in baited_digest_sets_chr.keys():
            n_total_interactions = 0
            n_di_interactions = 0
            n_uir_interactions = 0
            n_di_ne_interactions = 0
            n_uir_ne_interactions = 0
            for baited_digest in baited_digest_sets_chr[chr]:
                n_total_interactions += baited_digest.n_total_interactions()
                n_di_interactions += baited_digest.n_di_interactions()
                n_uir_interactions += baited_digest.n_uir_interactions()
                n_di_ne_interactions += baited_digest.n_di_ne_interactions()
                n_uir_ne_interactions += baited_digest.n_uir_ne_interactions()
            proportion_di = round(n_di_interactions/float(n_total_interactions), 2)
            proportion_uir = round(n_uir_interactions / float(n_total_interactions), 2)
            di_ne_proportion = round(n_di_ne_interactions/float(n_di_interactions), 2)
            uir_ne_proportion = round(n_uir_ne_interactions/float(n_uir_interactions), 2)
            print(chr + '\t' +
                  str(n_total_interactions) + '\t' +
                  str(n_di_interactions) + '\t' +
                  str(proportion_di) + '\t' +
                  str(di_ne_proportion) + '\t' +
                  str(proportion_uir) + '\t' +
                  str(uir_ne_proportion)
                  )
