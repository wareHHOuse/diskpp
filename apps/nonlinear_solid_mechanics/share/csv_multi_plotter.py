#!/usr/bin/env python3

"""Interactive multi-file CSV plotter for PyQt6 or PyQt5.

PyQt6 is attempted first. If it is unavailable, the application falls
back automatically to PyQt5.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

try:
    from PyQt6.QtCore import Qt
    from PyQt6.QtWidgets import (
        QApplication,
        QAbstractItemView,
        QCheckBox,
        QComboBox,
        QFileDialog,
        QGroupBox,
        QHBoxLayout,
        QLabel,
        QLineEdit,
        QListWidget,
        QListWidgetItem,
        QMainWindow,
        QMessageBox,
        QPushButton,
        QScrollArea,
        QSplitter,
        QSpinBox,
        QVBoxLayout,
        QWidget,
    )

    QT_API = "PyQt6"
    HORIZONTAL = Qt.Orientation.Horizontal
    TEXT_SELECTABLE_BY_MOUSE = Qt.TextInteractionFlag.TextSelectableByMouse
    ITEM_IS_USER_CHECKABLE = Qt.ItemFlag.ItemIsUserCheckable
    CHECKED = Qt.CheckState.Checked
    UNCHECKED = Qt.CheckState.Unchecked
    USER_ROLE = Qt.ItemDataRole.UserRole
    EXTENDED_SELECTION = QAbstractItemView.SelectionMode.ExtendedSelection
    MESSAGE_YES = QMessageBox.StandardButton.Yes
    MESSAGE_NO = QMessageBox.StandardButton.No
    CUSTOM_CONTEXT_MENU = Qt.ContextMenuPolicy.CustomContextMenu

except ImportError:
    from PyQt5.QtCore import Qt
    from PyQt5.QtWidgets import (
        QApplication,
        QAbstractItemView,
        QCheckBox,
        QComboBox,
        QFileDialog,
        QGroupBox,
        QHBoxLayout,
        QLabel,
        QLineEdit,
        QListWidget,
        QListWidgetItem,
        QMainWindow,
        QMessageBox,
        QPushButton,
        QScrollArea,
        QSplitter,
        QSpinBox,
        QVBoxLayout,
        QWidget,
    )

    QT_API = "PyQt5"
    HORIZONTAL = Qt.Horizontal
    TEXT_SELECTABLE_BY_MOUSE = Qt.TextSelectableByMouse
    ITEM_IS_USER_CHECKABLE = Qt.ItemIsUserCheckable
    CHECKED = Qt.Checked
    UNCHECKED = Qt.Unchecked
    USER_ROLE = Qt.UserRole
    EXTENDED_SELECTION = QAbstractItemView.ExtendedSelection
    MESSAGE_YES = QMessageBox.Yes
    MESSAGE_NO = QMessageBox.No
    CUSTOM_CONTEXT_MENU = Qt.CustomContextMenu


from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qtagg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure


class CSVPlotter(QMainWindow):
    def __init__(self, filenames=None):
        super().__init__()
        self.dataframes = {}
        # User-defined legend labels, indexed by absolute file path.
        # Labels remain available when a file is unchecked and rechecked.
        self.legend_labels = {}
        self.setWindowTitle("CSV Results Comparison")
        self.resize(1500, 900)
        self._build_interface()
        if filenames:
            self.load_csv_files(filenames)

    def _build_interface(self):
        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        main_layout = QVBoxLayout(central_widget)

        top_layout = QHBoxLayout()
        self.open_button = QPushButton("Open CSV files")
        self.open_button.clicked.connect(self.select_files)
        self.clear_button = QPushButton("Close all")
        self.clear_button.clicked.connect(self.clear_files)
        self.summary_label = QLabel("No file loaded")
        self.summary_label.setTextInteractionFlags(TEXT_SELECTABLE_BY_MOUSE)
        top_layout.addWidget(self.open_button)
        top_layout.addWidget(self.clear_button)
        top_layout.addWidget(self.summary_label, stretch=1)
        main_layout.addLayout(top_layout)

        splitter = QSplitter(HORIZONTAL)
        main_layout.addWidget(splitter, stretch=1)

        control_scroll = QScrollArea()
        control_scroll.setWidgetResizable(True)
        control_scroll.setMinimumWidth(320)
        control_widget = QWidget()
        control_layout = QVBoxLayout(control_widget)
        control_scroll.setWidget(control_widget)
        splitter.addWidget(control_scroll)

        files_group = QGroupBox("Files to display")
        files_layout = QVBoxLayout(files_group)
        self.file_list = QListWidget()
        self.file_list.setSelectionMode(EXTENDED_SELECTION)
        self.file_list.setToolTip(
            "Left-click a file to edit its legend label. "
            "Right-click it to remove it."
        )
        self.file_list.itemChanged.connect(self.update_plot)
        self.file_list.currentItemChanged.connect(self._on_current_file_changed)
        self.file_list.setContextMenuPolicy(CUSTOM_CONTEXT_MENU)
        self.file_list.customContextMenuRequested.connect(
            self._remove_file_by_right_click
        )
        files_layout.addWidget(self.file_list)

        legend_editor_layout = QHBoxLayout()
        self.legend_label_caption = QLabel("Legend label:")
        self.legend_label_edit = QLineEdit()
        self.legend_label_edit.setPlaceholderText("Select a file to edit its legend label")
        self.legend_label_edit.setEnabled(False)
        self.legend_label_edit.editingFinished.connect(self._store_current_legend_label)
        legend_editor_layout.addWidget(self.legend_label_caption)
        legend_editor_layout.addWidget(self.legend_label_edit, stretch=1)
        files_layout.addLayout(legend_editor_layout)

        files_buttons_layout = QHBoxLayout()
        self.enable_all_files_button = QPushButton("All")
        self.enable_all_files_button.clicked.connect(self.enable_all_files)
        self.disable_all_files_button = QPushButton("None")
        self.disable_all_files_button.clicked.connect(self.disable_all_files)
        self.remove_files_button = QPushButton("Remove")
        self.remove_files_button.clicked.connect(self.remove_selected_files)
        files_buttons_layout.addWidget(self.enable_all_files_button)
        files_buttons_layout.addWidget(self.disable_all_files_button)
        files_buttons_layout.addWidget(self.remove_files_button)
        files_layout.addLayout(files_buttons_layout)
        control_layout.addWidget(files_group)

        x_group = QGroupBox("Common X axis")
        x_layout = QVBoxLayout(x_group)
        self.x_combo = QComboBox()
        self.x_combo.currentTextChanged.connect(self.update_plot)
        x_layout.addWidget(self.x_combo)
        control_layout.addWidget(x_group)

        variables_group = QGroupBox("Subplots to display")
        variables_layout = QVBoxLayout(variables_group)
        self.variable_list = QListWidget()
        self.variable_list.itemChanged.connect(self.update_plot)
        variables_layout.addWidget(self.variable_list)

        variables_buttons_layout = QHBoxLayout()
        self.select_all_button = QPushButton("All")
        self.select_all_button.clicked.connect(self.select_all_variables)
        self.select_none_button = QPushButton("None")
        self.select_none_button.clicked.connect(self.select_no_variables)
        variables_buttons_layout.addWidget(self.select_all_button)
        variables_buttons_layout.addWidget(self.select_none_button)
        variables_layout.addLayout(variables_buttons_layout)
        control_layout.addWidget(variables_group, stretch=1)

        layout_group = QGroupBox("Plot layout")
        layout_controls = QHBoxLayout(layout_group)

        layout_controls.addWidget(QLabel("Rows:"))
        self.rows_spinbox = QSpinBox()
        self.rows_spinbox.setRange(0, 20)
        self.rows_spinbox.setSpecialValueText("Auto")
        self.rows_spinbox.setValue(0)
        self.rows_spinbox.setToolTip(
            "Auto uses as many rows as required. If an explicit grid is too "
            "small, extra rows are added so that no selected subplot is hidden."
        )
        self.rows_spinbox.valueChanged.connect(self.update_plot)
        layout_controls.addWidget(self.rows_spinbox)

        layout_controls.addWidget(QLabel("Columns:"))
        self.columns_spinbox = QSpinBox()
        self.columns_spinbox.setRange(1, 20)
        self.columns_spinbox.setValue(1)
        self.columns_spinbox.valueChanged.connect(self.update_plot)
        layout_controls.addWidget(self.columns_spinbox)

        control_layout.addWidget(layout_group)

        options_group = QGroupBox("Options")
        options_layout = QVBoxLayout(options_group)
        self.grid_checkbox = QCheckBox("Grid")
        self.grid_checkbox.setChecked(True)
        self.grid_checkbox.stateChanged.connect(self.update_plot)
        self.legend_checkbox = QCheckBox("Legend")
        self.legend_checkbox.setChecked(True)
        self.legend_checkbox.stateChanged.connect(self.update_plot)
        self.log_x_checkbox = QCheckBox("Logarithmic X scale")
        self.log_x_checkbox.stateChanged.connect(self.update_plot)
        self.log_y_checkbox = QCheckBox("Symmetric logarithmic Y scale")
        self.log_y_checkbox.stateChanged.connect(self.update_plot)
        self.markers_checkbox = QCheckBox("Show markers")
        self.markers_checkbox.setChecked(False)
        self.markers_checkbox.stateChanged.connect(self.update_plot)
        options_layout.addWidget(self.grid_checkbox)
        options_layout.addWidget(self.legend_checkbox)
        options_layout.addWidget(self.log_x_checkbox)
        options_layout.addWidget(self.log_y_checkbox)
        options_layout.addWidget(self.markers_checkbox)
        control_layout.addWidget(options_group)

        self.export_button = QPushButton("Export figure")
        self.export_button.clicked.connect(self.export_figure)
        control_layout.addWidget(self.export_button)

        plot_widget = QWidget()
        plot_layout = QVBoxLayout(plot_widget)
        self.figure = Figure()
        self.canvas = FigureCanvas(self.figure)
        self.toolbar = NavigationToolbar(self.canvas, self)
        plot_layout.addWidget(self.toolbar)
        plot_layout.addWidget(self.canvas)
        splitter.addWidget(plot_widget)
        splitter.setSizes([350, 1150])
        self._draw_empty_figure()
        self.statusBar().showMessage(f"Ready ({QT_API})")

    @staticmethod
    def _read_csv(filename):
        dataframe = pd.read_csv(filename, sep=";", skipinitialspace=True, engine="python")
        dataframe.columns = [str(column).strip() for column in dataframe.columns]
        dataframe = dataframe.dropna(axis=1, how="all")
        for column in dataframe.columns:
            if dataframe[column].dtype == object:
                dataframe[column] = dataframe[column].astype(str).str.strip()
            converted = pd.to_numeric(dataframe[column], errors="coerce")
            if converted.notna().any():
                dataframe[column] = converted
        return dataframe

    @staticmethod
    def _numeric_columns(dataframe):
        return [
            column for column in dataframe.columns
            if pd.api.types.is_numeric_dtype(dataframe[column])
        ]

    def select_files(self):
        filenames, _ = QFileDialog.getOpenFileNames(
            self,
            "Open result files",
            "",
            "CSV files (*.csv *.txt);;All files (*)",
        )
        if filenames:
            self.load_csv_files(filenames)

    def load_csv_files(self, filenames):
        errors = []
        added = 0
        selected_variables = set(self.selected_variables())
        current_x = self.x_combo.currentText()

        for filename in filenames:
            path = str(Path(filename).resolve())
            if path in self.dataframes:
                continue
            try:
                dataframe = self._read_csv(path)
            except Exception as error:
                errors.append(f"{Path(path).name} : {error}")
                continue
            if dataframe.empty:
                errors.append(f"{Path(path).name}: empty file")
                continue
            if len(self._numeric_columns(dataframe)) < 2:
                errors.append(f"{Path(path).name}: fewer than two numeric columns")
                continue
            self.dataframes[path] = dataframe
            # Default legend label is the file name. It can later be edited.
            self.legend_labels.setdefault(path, Path(path).name)
            added += 1

        self._rebuild_file_list()
        self._rebuild_column_lists(selected_variables, current_x)
        self._update_summary()
        self.update_plot()

        if errors:
            QMessageBox.warning(
                self,
                "Some files could not be loaded",
                "\n".join(errors),
            )
        self.statusBar().showMessage(
            f"{added} file(s) added, {len(self.dataframes)} file(s) loaded"
        )

    def _rebuild_file_list(self):
        checked_paths = set(self.enabled_file_paths())
        self.file_list.blockSignals(True)
        self.file_list.clear()
        for path in self.dataframes:
            item = QListWidgetItem(Path(path).name)
            item.setFlags(item.flags() | ITEM_IS_USER_CHECKABLE)
            item.setData(USER_ROLE, path)
            item.setToolTip(path)
            item.setCheckState(
                CHECKED if not checked_paths or path in checked_paths else UNCHECKED
            )
            self.file_list.addItem(item)
        self.file_list.blockSignals(False)

        if self.file_list.count() > 0:
            self.file_list.setCurrentRow(0)
        else:
            self._on_current_file_changed(None, None)

    def enabled_file_paths(self):
        return [
            self.file_list.item(index).data(Qt.UserRole)
            for index in range(self.file_list.count())
            if self.file_list.item(index).checkState() == CHECKED
        ]

    def enable_all_files(self):
        self.file_list.blockSignals(True)
        for index in range(self.file_list.count()):
            self.file_list.item(index).setCheckState(CHECKED)
        self.file_list.blockSignals(False)
        self.update_plot()

    def disable_all_files(self):
        self.file_list.blockSignals(True)
        for index in range(self.file_list.count()):
            self.file_list.item(index).setCheckState(UNCHECKED)
        self.file_list.blockSignals(False)
        self.update_plot()

    def _remove_file_by_right_click(self, position):
        """Remove the file under the mouse after a right-click confirmation."""
        item = self.file_list.itemAt(position)

        # A right-click in an empty area of the list does nothing.
        if item is None:
            return

        path = item.data(USER_ROLE)
        filename = Path(path).name

        answer = QMessageBox.question(
            self,
            "Remove file",
            f'Remove "{filename}" from the application?',
            MESSAGE_YES | MESSAGE_NO,
            MESSAGE_NO,
        )

        if answer != MESSAGE_YES:
            return

        selected_variables = set(self.selected_variables())
        current_x = self.x_combo.currentText()

        self.dataframes.pop(path, None)
        self.legend_labels.pop(path, None)

        self._rebuild_file_list()
        self._rebuild_column_lists(selected_variables, current_x)
        self._update_summary()
        self.update_plot()


    def remove_selected_files(self):
        selected_items = self.file_list.selectedItems()
        if not selected_items:
            return
        selected_variables = set(self.selected_variables())
        current_x = self.x_combo.currentText()
        for item in selected_items:
            path = item.data(USER_ROLE)
            self.dataframes.pop(path, None)
            self.legend_labels.pop(path, None)
        self._rebuild_file_list()
        self._rebuild_column_lists(selected_variables, current_x)
        self._update_summary()
        self.update_plot()

    def clear_files(self):
        self.dataframes.clear()
        self.legend_labels.clear()
        self.file_list.clear()
        self.legend_label_edit.clear()
        self.legend_label_edit.setEnabled(False)
        self.x_combo.clear()
        self.variable_list.clear()
        self._update_summary()
        self._draw_empty_figure()
        self.statusBar().showMessage("All files have been closed")

    def _all_numeric_columns(self):
        columns = set()
        for dataframe in self.dataframes.values():
            columns.update(self._numeric_columns(dataframe))
        return sorted(columns)

    def _common_numeric_columns(self):
        if not self.dataframes:
            return []
        column_sets = [
            set(self._numeric_columns(dataframe)) for dataframe in self.dataframes.values()
        ]
        return sorted(set.intersection(*column_sets))

    def _rebuild_column_lists(self, selected_variables=None, preferred_x=None):
        selected_variables = selected_variables or set()
        all_columns = self._all_numeric_columns()
        common_columns = self._common_numeric_columns()

        self.x_combo.blockSignals(True)
        self.x_combo.clear()
        self.x_combo.addItems(common_columns)
        if preferred_x in common_columns:
            self.x_combo.setCurrentText(preferred_x)
        elif "Time" in common_columns:
            self.x_combo.setCurrentText("Time")
        elif common_columns:
            self.x_combo.setCurrentIndex(0)
        self.x_combo.blockSignals(False)

        x_column = self.x_combo.currentText()
        default_variables = {"ux", "uy", "vx", "vy", "ax", "ay", "sxx", "syy", "szz", "sxy"}

        self.variable_list.blockSignals(True)
        self.variable_list.clear()
        for column in all_columns:
            if column == x_column:
                continue
            item = QListWidgetItem(column)
            item.setFlags(item.flags() | ITEM_IS_USER_CHECKABLE)
            checked = column in selected_variables if selected_variables else column in default_variables
            item.setCheckState(CHECKED if checked else UNCHECKED)
            if column not in common_columns:
                existing_count = sum(column in df.columns for df in self.dataframes.values())
                item.setToolTip(
                    f"Column available in {existing_count}/{len(self.dataframes)} file(s)"
                )
            self.variable_list.addItem(item)
        self.variable_list.blockSignals(False)

    def selected_variables(self):
        return [
            self.variable_list.item(index).text()
            for index in range(self.variable_list.count())
            if self.variable_list.item(index).checkState() == CHECKED
        ]

    def select_all_variables(self):
        self.variable_list.blockSignals(True)
        for index in range(self.variable_list.count()):
            self.variable_list.item(index).setCheckState(CHECKED)
        self.variable_list.blockSignals(False)
        self.update_plot()

    def select_no_variables(self):
        self.variable_list.blockSignals(True)
        for index in range(self.variable_list.count()):
            self.variable_list.item(index).setCheckState(UNCHECKED)
        self.variable_list.blockSignals(False)
        self.update_plot()

    def _draw_empty_figure(self, message=None):
        self.figure.clear()
        axis = self.figure.add_subplot(111)
        axis.set_axis_off()
        axis.text(
            0.5,
            0.5,
            message or "Open one or more CSV files",
            ha="center",
            va="center",
            transform=axis.transAxes,
            fontsize=13,
        )
        self.canvas.draw_idle()

    def _on_current_file_changed(self, current, previous):
        """Load the selected file's persistent legend label into the editor."""
        if current is None:
            self.legend_label_edit.clear()
            self.legend_label_edit.setEnabled(False)
            return

        path = current.data(Qt.UserRole)
        self.legend_label_edit.blockSignals(True)
        self.legend_label_edit.setText(
            self.legend_labels.get(path, Path(path).name)
        )
        self.legend_label_edit.blockSignals(False)
        self.legend_label_edit.setEnabled(True)

    def _store_current_legend_label(self):
        """Store the edited label without depending on the checked state."""
        item = self.file_list.currentItem()
        if item is None:
            return

        path = item.data(USER_ROLE)
        label = self.legend_label_edit.text().strip()

        if not label:
            label = Path(path).name
            self.legend_label_edit.setText(label)

        self.legend_labels[path] = label
        self.update_plot()

    def _plot_labels(self, paths):
        """Return persistent user-defined labels for the requested paths."""
        return {
            path: self.legend_labels.get(path, Path(path).name)
            for path in paths
        }

    def update_plot(self):
        if not self.dataframes:
            self._draw_empty_figure()
            return
        enabled_paths = self.enabled_file_paths()
        if not enabled_paths:
            self._draw_empty_figure("Enable at least one file")
            return
        x_column = self.x_combo.currentText()
        if not x_column:
            self._draw_empty_figure("No common X axis is available")
            return
        selected_variables = [
            variable for variable in self.selected_variables() if variable != x_column
        ]
        if not selected_variables:
            self._draw_empty_figure("Select at least one variable")
            return

        number_of_plots = len(selected_variables)
        requested_rows = self.rows_spinbox.value()
        number_of_columns = self.columns_spinbox.value()

        # Rows == 0 means automatic. The default (Auto, 1 column) therefore
        # preserves the original vertical layout.
        required_rows = (number_of_plots + number_of_columns - 1) // number_of_columns
        if requested_rows == 0:
            number_of_rows = required_rows
        else:
            # Never silently omit a selected variable when the requested grid
            # is too small. Add the minimum number of rows required.
            number_of_rows = max(requested_rows, required_rows)

        self.figure.clear()
        axes_grid = self.figure.subplots(
            nrows=number_of_rows,
            ncols=number_of_columns,
            sharex=True,
            squeeze=False,
        )
        axes = axes_grid.ravel()
        active_axes = axes[:number_of_plots]
        empty_axes = axes[number_of_plots:]

        # Keep unused grid cells completely blank.
        for axis in empty_axes:
            axis.set_axis_off()

        labels = self._plot_labels(enabled_paths)
        marker = "o" if self.markers_checkbox.isChecked() else None
        plotted_curve_count = 0

        for axis, variable in zip(active_axes, selected_variables):
            variable_curve_count = 0
            for path in enabled_paths:
                dataframe = self.dataframes[path]
                if x_column not in dataframe.columns or variable not in dataframe.columns:
                    continue
                x_values = dataframe[x_column].to_numpy(dtype=float)
                y_values = dataframe[variable].to_numpy(dtype=float)
                valid = np.isfinite(x_values) & np.isfinite(y_values)
                if not np.any(valid):
                    continue
                axis.plot(
                    x_values[valid],
                    y_values[valid],
                    linewidth=1.5,
                    marker=marker,
                    markersize=3.0,
                    markevery=max(1, np.count_nonzero(valid) // 30),
                    label=labels[path],
                )
                variable_curve_count += 1
                plotted_curve_count += 1

            axis.set_ylabel(variable)
            axis.grid(self.grid_checkbox.isChecked(), which="both", alpha=0.35)
            axis.set_xscale("log" if self.log_x_checkbox.isChecked() else "linear")
            if self.log_y_checkbox.isChecked():
                axis.set_yscale("symlog", linthresh=1.0e-12)
            else:
                axis.set_yscale("linear")
            if self.legend_checkbox.isChecked() and variable_curve_count > 0:
                axis.legend(loc="best", fontsize="small")
            if variable_curve_count == 0:
                axis.text(
                    0.5,
                    0.5,
                    "Variable missing from enabled files",
                    ha="center",
                    va="center",
                    transform=axis.transAxes,
                )

        # Each populated subplot gets an X label because, with multiple
        # columns, the final populated cell may not be in the last grid row.
        for axis in active_axes:
            axis.set_xlabel(x_column)

        active_axes[0].set_title("Results comparison")
        self.figure.set_size_inches(
            max(10.0, 5.0 * number_of_columns),
            max(6.0, 2.8 * number_of_rows),
            forward=True,
        )
        self.figure.tight_layout()
        self.canvas.draw_idle()
        self.statusBar().showMessage(
            f"{len(enabled_paths)} active file(s), "
            f"{len(selected_variables)} subplot(s), "
            f"grid {number_of_rows} x {number_of_columns}, "
            f"{plotted_curve_count} curve(s)"
        )

    def export_figure(self):
        if not self.dataframes:
            QMessageBox.information(self, "No data", "Load at least one CSV file first.")
            return
        filename, selected_filter = QFileDialog.getSaveFileName(
            self,
            "Export figure",
            "comparaison.png",
            "PNG image (*.png);;PDF document (*.pdf);;SVG image (*.svg)",
        )
        if not filename:
            return
        if not Path(filename).suffix:
            if "PDF" in selected_filter:
                filename += ".pdf"
            elif "SVG" in selected_filter:
                filename += ".svg"
            else:
                filename += ".png"
        try:
            self.figure.savefig(filename, dpi=300, bbox_inches="tight")
        except Exception as error:
            QMessageBox.critical(
                self, "Export error", f"Could not export the figure:\n\n{error}"
            )
            return
        self.statusBar().showMessage(f"Figure exported: {filename}")

    def _update_summary(self):
        number_of_files = len(self.dataframes)
        if number_of_files == 0:
            self.summary_label.setText("No file loaded")
            return
        total_rows = sum(len(dataframe) for dataframe in self.dataframes.values())
        self.summary_label.setText(
            f"{number_of_files} file(s), {total_rows} total row(s)"
        )


def main():
    application = QApplication(sys.argv)
    application.setApplicationName("CSV Multi-File Plotter")
    window = CSVPlotter(sys.argv[1:])
    window.show()
    if QT_API == "PyQt6":
        sys.exit(application.exec())
    sys.exit(application.exec_())


if __name__ == "__main__":
    main()
