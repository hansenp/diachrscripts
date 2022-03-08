from matplotlib.ticker import FuncFormatter
import matplotlib.pyplot as plt
from scipy import stats


def thousands(x, pos):
    return '%1.0fk' % (x*1e-3)


def millions(x, pos):
    return '%1.1fM' % (x*1e-6)


class IaFreqDistAnalysis_2:

    def __init__(self, input_table_path: str = None):

        # Dictionary that contains the numbers grouped by cell type, interaction and enrichment categories
        self._grouped_numbers = dict()

        self._names = {
            'DIX': 'Unbalanced without reference',
            'DI': 'Unbalanced',
            'UIR': 'Balanced reference',
            'UI': 'Balanced',
            'ALL': 'All',
            'N': 'Interaction numbers',
            'MED': 'Median interaction distances',
            'MAD': 'Median absolute deviations'
        }

        self._colors = {
            'DIX': 'orangered',
            'DI': 'orange',
            'UIR': 'lightblue',
            'UI': 'lightgray',
            'ALL': 'cornflowerblue'
        }

        # Read table
        with open(input_table_path, 'rt') as fp:
            header_row_dict = dict()
            for line in fp:

                fields = line.rstrip('\n').split('\t')
                if len(fields) != 63:
                    print('[ERROR] Malformed input line:\n\t' + line)
                    exit(1)

                if fields[0] == 'DESCRIPTION':  # Header row
                    self._grouped_numbers['CELL_TYPES'] = []
                    col_num = 0
                    for field in fields:
                        header_row_dict[col_num] = field
                        if 1 < col_num < 62:  # Init lists of numbers
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
                        if 1 < col_num < 62:  # Append numbers
                            i_cat, e_cat, m_type = header_row_dict[col_num].split('|')
                            self._grouped_numbers[i_cat][e_cat][m_type].append(int(field))
                        col_num += 1

    def create_i_dist_scatterplot(self,
                           i_cat_1: str = None,
                           i_cat_2: str = None,
                           e_cat_1: str = None,
                           e_cat_2: str = None,
                           m_type: str = None,
                           pdf_file_name: str = 'scatter_plot.pdf'):

        # Collect plot data
        cell_types = self._grouped_numbers['CELL_TYPES']
        x = self._grouped_numbers[i_cat_1][e_cat_1][m_type]
        y = self._grouped_numbers[i_cat_2][e_cat_2][m_type]
        m_type_name = self._names[m_type]
        i_cat_1_name = self._names[i_cat_1]
        i_cat_2_name = self._names[i_cat_2]
        i_cat_1_color = self._colors[i_cat_1]
        i_cat_2_color = self._colors[i_cat_2]

        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(4, 4))

        x_min = min(x)
        x_max = max(x)
        y_min = min(y)
        y_max = max(y)
        xy_min = min(x_min, y_min)
        xy_max = max(x_max, y_max)
        ax.scatter(x,
                   y,
                   c=[i_cat_1_color] * 17,
                   edgecolors=[i_cat_2_color] * 17,
                   linewidth=2)

        # Add title and axis labels
        ax.set_title(m_type_name + ' [bp]', loc='left')
        ax.set_xlabel(i_cat_1_name + ' - ' + e_cat_1)
        ax.set_ylabel(i_cat_2_name + ' - ' + e_cat_2)

        # Same number of ticks
        ax.locator_params(axis="x", nbins=4)
        ax.locator_params(axis="y", nbins=4)

        # Format ticks - thousands or millions
        if xy_max < 1000000:
            formatter = FuncFormatter(thousands)
        else:
            formatter = FuncFormatter(millions)
        ax.xaxis.set_major_formatter(formatter)
        ax.yaxis.set_major_formatter(formatter)

        # Rotate labels by 45 degrees
        ax.tick_params(labelrotation=45)

        # Add diagonal line
        ax.plot([xy_min, xy_max], [xy_min, xy_max], color='black', linewidth=0.5, linestyle='-.')

        # Add cell type names to scatterplot
        for i, cell_type in enumerate(cell_types):
            ax.annotate(cell_type, (x[i], y[i]), size=4)

        # Calculate P-value using Wilcoxon signed-rank test
        w_stat, p_value = stats.wilcoxon(x=x, y=y)

        # Calculate average deviation from zero
        diffs = [(x[i] - y[i]) for i in range(0, len(x))]
        avg_dev_zero = sum(diffs) / len(diffs)

        # Print test results
        print('P-value: ' + "{:.2e}".format(p_value))
        print('Test statistic: ' + "{:.2f}".format(w_stat))
        print('Average difference of medians: ' + "{:,.0f}".format(avg_dev_zero))

        # Save and return figure
        ax.set_aspect(1)
        fig.tight_layout()
        fig.savefig(pdf_file_name)
        return fig
