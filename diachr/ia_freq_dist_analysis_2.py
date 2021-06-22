import matplotlib.pyplot as plt

class IaFreqDistAnalysis_2:

    def __init__(self, input_table_path: str = None):

        # Dictionary that contains the numbers grouped by cell type, interaction and enrichment categories
        self._grouped_numbers = dict()

        self._names = {
            'DI': 'Directed',
            'DI_S': 'Directed simple',
            'DI_T': 'Directed twisted',
            'UIR': 'Undirected reference',
            'UI': 'Undirected',
            'ALL': 'All',
            'N': 'Interaction number',
            'MED': 'Median interaction distance',
            'MAD': 'Median absolute deviation'
        }

        self._colors = {
            'DI': 'orange',
            'DI_S': 'pink',
            'DI_T': 'cadetblue',
            'UIR': 'lightblue',
            'UI': 'lightgray',
            'ALL': 'cornflowerblue'
        }

        # Read table
        with open(input_table_path, 'rt') as fp:
            header_row_dict = dict()
            for line in fp:

                fields = line.rstrip('\n').split('\t')
                if len(fields) != 75:
                    print('[ERROR] Malformed input line:\n\t' + line)
                    exit(1)

                if fields[0] == 'DESCRIPTION': # Header row
                    self._grouped_numbers['CELL_TYPES'] = []
                    col_num = 0
                    for field in fields:
                        header_row_dict[col_num] = field
                        if 1 < col_num < 74: # Init lists of numbers
                            i_cat, e_cat, m_type = field.split('|')
                            if i_cat not in self._grouped_numbers:
                                self._grouped_numbers[i_cat] = dict()
                            if e_cat not in self._grouped_numbers[i_cat]:
                                self._grouped_numbers[i_cat][e_cat] = dict()
                            if m_type not in self._grouped_numbers[i_cat][e_cat]:
                                self._grouped_numbers[i_cat][e_cat][m_type] = []
                        col_num += 1
                else:
                    col_num = 0
                    for field in fields:
                        if col_num == 1:
                            self._grouped_numbers['CELL_TYPES'].append(field)
                        if 1 < col_num < 74:  # Append numbers
                            i_cat, e_cat, m_type = header_row_dict[col_num].split('|')
                            self._grouped_numbers[i_cat][e_cat][m_type].append(int(field))
                        col_num += 1

    def create_scatterplot(self,
                           i_cat_1: str = None,
                           i_cat_2: str = None,
                           e_cat_1: str = None,
                           e_cat_2: str = None,
                           m_type: str = None,
                           pdf_file_name: str = 'scatter_plot.pdf'):


        x = self._grouped_numbers[i_cat_1][e_cat_1][m_type]
        y = self._grouped_numbers[i_cat_2][e_cat_2][m_type]

        m_type = self._names[m_type]
        i_cat_1_name = self._names[i_cat_1]
        i_cat_2_name = self._names[i_cat_2]
        i_cat_1_color = self._colors[i_cat_1]
        i_cat_2_color = self._colors[i_cat_2]

        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(4, 4))

        cell_types = self._grouped_numbers['CELL_TYPES']

        x_min = min(x)
        x_max = max(x)
        y_min = min(y)
        y_max = max(y)
        ax.scatter(x,
                   y,
                   c = [i_cat_1_color]*17,
                   edgecolors = [i_cat_2_color]*17)

        ax.set_title(m_type)
        ax.set_xlabel(i_cat_1_name + ' - ' + e_cat_1)
        ax.set_ylabel(i_cat_2_name + ' - ' + e_cat_2)
        xy_ticks, xy_tick_labels = self.make_ticks(max(x_max, y_max))

        ax.set_xticks(xy_ticks)
        ax.set_xticklabels(xy_tick_labels)
        ax.set_yticks(xy_ticks)
        ax.set_yticklabels(xy_tick_labels)
        ax.plot([min(x_min, y_min), max(x_max, y_max)], [min(x_min, y_min), max(x_max, y_max)])

        # Save and return figure
        fig.tight_layout(pad=1.25)
        fig.savefig(pdf_file_name)
        return fig

    def make_ticks(self, max_val: int = 100):

        if max_val < 10:
            return list(range(0, 50 + 1, 1)), list(range(0, 50 + 1, 1))

        if max_val < 25:
            return list(range(0, 50 + 1, 5)), list(range(0, 50 + 1, 5))

        if max_val < 50:
            return list(range(0, 50 + 1, 10)), list(range(0, 50 + 1, 10))

        # Create list of maximal ticks
        tick_lims = []
        x = [10, 25, 50]
        while x[2] <= max_val:
            x = [(y * 10) for y in x]
            tick_lims = tick_lims + x

        # Find maximal tick for input value
        tick_max_idx = 0
        while tick_lims[tick_max_idx] <= max_val:
            tick_max_idx += 1

        ticks = list(range(0, int(tick_lims[tick_max_idx]) + 1, int(tick_lims[tick_max_idx] / 5)))

        if max_val < 1000:
            tick_labels = ticks
        elif max_val < 1000000:
            tick_labels = []
            for tick in ticks:
                tick_labels.append(str(tick / 1000) + 'k')
        else:
            tick_labels = []
            for tick in ticks:
                tick_labels.append(str(tick / 1000000) + 'M')

        return ticks, tick_labels
