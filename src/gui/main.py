import gi
gi.require_version('Gtk', '3.0')
from gi.repository import Gtk
import json
import re
import sys

MAX_VAL = 999999

class App(Gtk.Window):
    def __init__(self):
        super().__init__()
        self.set_default_size(950, 750)
        self.connect("destroy", Gtk.main_quit)

        self.grid_list = []
        self.terms_data = []

        main_box = Gtk.Box(orientation=Gtk.Orientation.VERTICAL)
        self.add(main_box)

        #region MENU_BAR
        mb = Gtk.MenuBar()

        filemenu = Gtk.Menu()
        filem = Gtk.MenuItem(label="File")
        filem.set_submenu(filemenu)

        self.import_btn = Gtk.MenuItem(label="Import settings...")
        self.import_btn.connect("activate", self.on_import)
        filemenu.append(self.import_btn)

        self.export_btn = Gtk.MenuItem(label="Export settings...")
        #self.export_btn.get_style_context().add_class("suggested-action")
        self.export_btn.connect("activate", self.on_export)
        filemenu.append(self.export_btn)

        mb.append(filem)
        main_box.add(mb)
        #endregion

        body_box = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=10)
        body_box.set_margin_start(10)
        body_box.set_margin_end(10)
        body_box.set_margin_top(10)
        body_box.set_margin_bottom(10)
        main_box.add(body_box)

        #region SETTINGS
        ctrl_panel = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=10)
        body_box.add(ctrl_panel)

        settings_frame = Gtk.Frame(label="Term Library Settings")
        settings_box = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=8)
        settings_box.set_margin_start(10)
        settings_box.set_margin_end(10)
        settings_box.set_margin_top(10)
        settings_box.set_margin_bottom(10)
        settings_frame.add(settings_box)
        ctrl_panel.pack_start(settings_frame, False, False, 0)

        grid_g = Gtk.Grid()
        grid_g.set_column_spacing(6)
        grid_g.set_row_spacing(6)
        settings_box.pack_start(grid_g, False, False, 0)

        grid_g.attach(Gtk.Label(label="Dimensions"), 0, 0, 1, 1)
        self.grid_entry = Gtk.Entry()
        self.grid_entry.connect("changed", self.on_settings_changed)
        grid_g.attach(self.grid_entry, 1, 0, 1, 1)

        grid_g.attach(Gtk.Label(label="Fields"), 0, 1, 1, 1)
        self.fields_entry = Gtk.Entry()
        self.fields_entry.connect("changed", self.on_settings_changed)
        grid_g.attach(self.fields_entry, 1, 1, 1, 1)

        grid_g.attach(Gtk.Label(label="Size of Data"), 0, 2, 1, 1)
        self.size_data_entry = Gtk.Entry()
        self.size_data_entry.connect("changed", self.on_settings_changed)
        grid_g.attach(self.size_data_entry, 1, 2, 1, 1)

        grid_g.attach(Gtk.Label(label="Integration Size"), 0, 3, 1, 1)
        self.int_size_entry = Gtk.Entry()
        self.int_size_entry.connect("changed", self.on_settings_changed)
        grid_g.attach(self.int_size_entry, 1, 3, 1, 1)

        grid_g.attach(Gtk.Label(label="Seed"), 0, 4, 1, 1)
        self.seed_spin = Gtk.SpinButton.new_with_range(0, MAX_VAL, 1)
        grid_g.attach(self.seed_spin, 1, 4, 1, 1)

        grid_g.attach(Gtk.Label(label="Gamma"), 0, 5, 1, 1)
        self.gamma_spin = Gtk.SpinButton.new_with_range(0.0, float(MAX_VAL), 0.1)
        self.gamma_spin.set_digits(8)
        grid_g.attach(self.gamma_spin, 1, 5, 1, 1)

        grid_g.attach(Gtk.Label(label="# Equations"), 0, 6, 1, 1)
        self.eq_spin = Gtk.SpinButton.new_with_range(1, MAX_VAL, 1)
        grid_g.attach(self.eq_spin, 1, 6, 1, 1)

        grid_g.attach(Gtk.Label(label="# Windows"), 0, 7, 1, 1)
        self.windows_spin = Gtk.SpinButton.new_with_range(1, MAX_VAL, 1)
        grid_g.attach(self.windows_spin, 1, 7, 1, 1)

        self.on_settings_changed(None)
        #endregion

        #region BULK_GENERATOR
        bulk_frame = Gtk.Frame(label="Bulk Derivative Generator")
        bulk_box = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=8)
        bulk_box.set_margin_start(10)
        bulk_box.set_margin_end(10)
        bulk_box.set_margin_top(10)
        bulk_box.set_margin_bottom(10)
        bulk_frame.add(bulk_box)
        ctrl_panel.pack_start(bulk_frame, False, False, 0)

        bulk_grid = Gtk.Grid()
        bulk_grid.set_column_spacing(6)
        bulk_grid.set_row_spacing(6)
        bulk_box.pack_start(bulk_grid, False, False, 0)

        bulk_grid.attach(Gtk.Label(label="Base Term"), 0, 0, 1, 1)
        self.bulk_base_entry = Gtk.Entry()
        self.bulk_base_entry.set_text("u")
        bulk_grid.attach(self.bulk_base_entry, 1, 0, 1, 1)

        bulk_grid.attach(Gtk.Label(label="Order"), 0, 1, 1, 1)
        self.bulk_k_spin = Gtk.SpinButton.new_with_range(0, 99, 1)
        bulk_grid.attach(self.bulk_k_spin, 1, 1, 1, 1)

        bulk_grid.attach(Gtk.Label(label="Scale"), 0, 2, 1, 1)
        self.bulk_scale_spin = Gtk.SpinButton.new_with_range(-float(MAX_VAL), float(MAX_VAL), 0.01)
        self.bulk_scale_spin.set_digits(8)
        self.bulk_scale_spin.set_value(1.0)
        bulk_grid.attach(self.bulk_scale_spin, 1, 2, 1, 1)

        self.bulk_btn = Gtk.Button(label="Generate")
        self.bulk_btn.get_style_context().add_class("suggested-action")
        self.bulk_btn.connect("clicked", self.on_generate_bulk)
        bulk_grid.attach(self.bulk_btn, 0, 3, 4, 1)
        #endregion
        
        #region TERMS
        terms_frame = Gtk.Frame()
        terms_vbox = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=8)
        terms_vbox.set_margin_start(10)
        terms_vbox.set_margin_end(10)
        terms_vbox.set_margin_top(10)
        terms_vbox.set_margin_bottom(10)
        terms_frame.add(terms_vbox)
        body_box.pack_start(terms_frame, True, True, 0)

        input_hbox = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=6)
        terms_vbox.pack_start(input_hbox, False, False, 0)

        self.key_entry = Gtk.Entry()
        self.key_entry.set_hexpand(True)
        self.key_entry.set_placeholder_text("Term...")
        self.key_entry.connect("activate", lambda w: self.add_term())
        input_hbox.pack_start(self.key_entry, True, True, 0)

        self.scale_spin = Gtk.SpinButton.new_with_range(0.0, float(MAX_VAL), 0.1)
        self.scale_spin.set_digits(5)
        self.scale_spin.set_placeholder_text("Scale...")
        input_hbox.pack_start(self.scale_spin, False, False, 0)

        self.partial_btn = Gtk.Button(label="∂")
        self.partial_btn.connect("clicked", self.on_partial_clicked)
        input_hbox.pack_start(self.partial_btn, False, False, 0)

        self.add_btn = Gtk.Button(label="Add Term")
        self.add_btn.get_style_context().add_class("suggested-action")
        self.add_btn.connect("clicked", lambda w: self.add_term())
        input_hbox.pack_start(self.add_btn, False, False, 0)

        self.tree_store = Gtk.ListStore(str, str, str)
        self.tree_view = Gtk.TreeView(model=self.tree_store)
        self.tree_view.set_headers_visible(True)
        self.tree_view.set_grid_lines(Gtk.TreeViewGridLines.BOTH)
        for i, title in enumerate(["Key", "Scale", "Raw"]):
            renderer = Gtk.CellRendererText()
            column = Gtk.TreeViewColumn(title, renderer, text=i)
            if i == 2: column.set_expand(True)
            self.tree_view.append_column(column)

        scroll = Gtk.ScrolledWindow()
        scroll.set_policy(Gtk.PolicyType.AUTOMATIC, Gtk.PolicyType.AUTOMATIC)
        scroll.add(self.tree_view)
        terms_vbox.pack_start(scroll, True, True, 0)

        ctrl_hbox = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=6)
        terms_vbox.pack_start(ctrl_hbox, False, False, 0)
        
        self.remove_btn = Gtk.Button(label="Remove Selected")
        self.remove_btn.connect("clicked", self.on_remove_term)
        ctrl_hbox.pack_start(self.remove_btn, False, False, 0)
        
        self.clear_btn = Gtk.Button(label="Clear All")
        self.clear_btn.connect("clicked", self.on_clear_terms)
        ctrl_hbox.pack_start(self.clear_btn, False, False, 0)
        #endregion

        #region MODEL_DISCOVERY

        #endregion

    def on_partial_clicked(self, button):
        pos = self.key_entry.get_position()
        text = self.key_entry.get_text()
        self.key_entry.set_text(text[:pos] + '∂' + text[pos:])
        self.key_entry.set_position(pos + 1)
        self.key_entry.grab_focus()

    def on_settings_changed(self, widget):
        self.grid_list = [x.strip() for x in self.grid_entry.get_text().split(',') if x.strip()]

    def parse_numerator(self, num):
        num = num.strip()
        if num.startswith('∂'):
            num = num[1:].strip()

        num = re.sub(r'^\^(\d+)\s*', '', num).strip()
        if not num: return None

        m = re.match(r'^(.+?)\^2$', num)
        if m: return [{".^2": [m.group(1).strip()]}]

        m = re.match(r'^(sin|cos|tan|exp)[.\s]*\((.+)\)$', num, re.IGNORECASE)
        if m: return [{f"{m.group(1).lower()}.": [m.group(2).strip()]}]

        for op_sym, json_op in [('+', '.+'), ('-', '.-'), ('*', '.*')]:
            if op_sym in num:
                parts = num.split(op_sym, 1)
                if len(parts) == 2 and parts[0].strip() and parts[1].strip():
                    return [{json_op: [p.strip() for p in parts]}]

        if '^' in num:
            parts = num.split('^', 1)
            if len(parts) == 2 and parts[0].strip() and parts[1].strip() and parts[1].strip() != '2':
                return [{".^n": [p.strip() for p in parts]}]

        return [{" ": [num]}]

    def parse_term_string(self, key, grid):
        key = key.strip()
        if '/' not in key: return None
        num_raw, den_raw = key.split('/', 1)

        den_match = re.match(r'∂([a-zA-Z0-9_]+)(?:\^(\d+))?', den_raw.strip())
        if not den_match: return None
        grid_var = den_match.group(1)
        deriv_order = int(den_match.group(2)) if den_match.group(2) else 1
        try:
            idx = grid.index(grid_var) + 1
        except ValueError:
            idx = 1
        derivs = [idx] * deriv_order

        term = self.parse_numerator(num_raw)
        if not term: return None

        return {"term": term, "derivs": derivs}

    def add_term(self):
        key = self.key_entry.get_text().strip()
        scale = round(self.scale_spin.get_value(), 4)
        if not key: return

        parsed = self.parse_term_string(key, self.grid_list)
        if not parsed:
            dialog = Gtk.MessageDialog(transient_for=self, flags=0,
                message_type=Gtk.MessageType.WARNING, buttons=Gtk.ButtonsType.OK,
                text="Invalid term")
            dialog.run()
            dialog.destroy()
            return

        self.terms_data.append({"key": key, "scale": scale, "term": parsed["term"], "derivs": parsed["derivs"]})
        preview = json.dumps({"term": parsed["term"], "derivs": parsed["derivs"]}, separators=(',', ':'))
        self.tree_store.append([key, str(scale), preview])
        self.key_entry.set_text("")
        self.key_entry.grab_focus()

    def on_generate_bulk(self, button):
        base = self.bulk_base_entry.get_text().strip()
        order = int(self.bulk_k_spin.get_value())
        scale = round(self.bulk_scale_spin.get_value(), 4)
        grid = self.grid_list

        if not base or not grid: return

        added_count = 0
        for g in grid:
            order_prefix = f"^{order}" if order > 1 else ""
            den_power = f"^{order}" if order > 1 else ""
            key = f"∂{order_prefix}{base}/∂{g}{den_power}"
            
            if any(t["key"] == key for t in self.terms_data):
                continue

            parsed = self.parse_term_string(key, grid)
            if parsed:
                self.terms_data.append({"key": key, "scale": scale, "term": parsed["term"], "derivs": parsed["derivs"]})
                preview = json.dumps({"term": parsed["term"], "derivs": parsed["derivs"]}, separators=(',', ':'))
                self.tree_store.append([key, str(scale), preview])
                added_count += 1

        dialog = Gtk.MessageDialog(transient_for=self, flags=0,
            message_type=Gtk.MessageType.INFO, buttons=Gtk.ButtonsType.OK,
            text=f"Added {added_count} terms")
        dialog.run()
        dialog.destroy()

    def on_remove_term(self, button):
        selection = self.tree_view.get_selection()
        model, treeiter = selection.get_selected()
        if treeiter:
            path = model.get_path(treeiter)
            idx = path.get_indices()[0]
            model.remove(treeiter)
            if 0 <= idx < len(self.terms_data):
                self.terms_data.pop(idx)

    def on_clear_terms(self, button):
        self.tree_store.clear()
        self.terms_data.clear()

    def on_import(self, button):
        dialog = Gtk.FileChooserDialog(title="Import Configuration JSON", parent=self,
            action=Gtk.FileChooserAction.OPEN)
        dialog.add_buttons(Gtk.STOCK_CANCEL, Gtk.ResponseType.CANCEL, Gtk.STOCK_OPEN, Gtk.ResponseType.OK)
        f_filter = Gtk.FileFilter()
        f_filter.set_name(".json")
        f_filter.add_mime_type("application/json")
        dialog.add_filter(f_filter)

        response = dialog.run()
        if response == Gtk.ResponseType.OK:
            filepath = dialog.get_filename()
            try:
                with open(filepath, 'r', encoding='utf-8') as f:
                    data = json.load(f)

                self.grid_entry.set_text(", ".join(data.get("grid", [])))
                fields = [k.strip() for k in data.get("fields", {}).keys()]
                self.fields_entry.set_text(", ".join(fields))
                self.size_data_entry.set_text(", ".join(map(str, data.get("size_of_data", []))))
                self.int_size_entry.set_text(", ".join(map(str, data.get("integration_size", []))))
                self.seed_spin.set_value(data.get("seed", 42))
                self.gamma_spin.set_value(data.get("gamma", 3.5))
                self.eq_spin.set_value(data.get("equations_to_discover", 5))
                self.windows_spin.set_value(data.get("number_of_windows", 64))
                self.on_settings_changed(None)

                self.tree_store.clear()
                self.terms_data.clear()
                lib_terms = data.get("library_terms", {})
                for key, val in lib_terms.items():
                    key = key.strip()
                    term = val.get("term", [])
                    derivs = val.get("derivs", [])
                    scale = val.get("scale", 1.0)
                    self.terms_data.append({"key": key, "scale": scale, "term": term, "derivs": derivs})
                    preview = json.dumps({"term": term, "derivs": derivs}, separators=(',', ':'))
                    self.tree_store.append([key, str(scale), preview])
            except Exception as e:
                err_dialog = Gtk.MessageDialog(transient_for=self, flags=0,
                    message_type=Gtk.MessageType.ERROR, buttons=Gtk.ButtonsType.OK,
                    text="Import Failed")
                err_dialog.format_secondary_text(str(e))
                err_dialog.run()
                err_dialog.destroy()
        dialog.destroy()

    def on_export(self, button):
        fields_dict = {f.strip(): f"result/value/{f.strip()}/value" for f in self.fields_entry.get_text().split(',') if f.strip()}
        try:
            size_data = [int(x.strip()) for x in self.size_data_entry.get_text().split(',') if x.strip()]
            int_size = [int(x.strip()) for x in self.int_size_entry.get_text().split(',') if x.strip()]
        except ValueError:
            size_data, int_size = [0], [0]

        library_terms = {}
        for td in self.terms_data:
            library_terms[td["key"]] = {
                "term": td["term"],
                "derivs": td["derivs"],
                "scale": td["scale"]
            }

        output = {
            "fields": fields_dict,
            "number_of_windows": int(self.windows_spin.get_value()),
            "degrees_of_freedom": 1,
            "dimensions": len(self.grid_list),
            "envelope_power": 4,
            "size_of_data": size_data,
            "integration_size": int_size,
            "buffer": 0,
            "seed": int(self.seed_spin.get_value()),
            "grid": self.grid_list,
            "gamma": self.gamma_spin.get_value(),
            "equations_to_discover": int(self.eq_spin.get_value()),
            "library_terms": library_terms
        }

        dialog = Gtk.FileChooserDialog(parent=self,
            action=Gtk.FileChooserAction.SAVE)
        dialog.add_buttons(Gtk.STOCK_CANCEL, Gtk.ResponseType.CANCEL, Gtk.STOCK_SAVE, Gtk.ResponseType.OK)
        dialog.set_current_name("settings.json")
        f_filter = Gtk.FileFilter()
        f_filter.set_name(".json")
        f_filter.add_mime_type("application/json")
        dialog.add_filter(f_filter)

        response = dialog.run()
        if response == Gtk.ResponseType.OK:
            filepath = dialog.get_filename()
            with open(filepath, 'w', encoding='utf-8') as f:
                json.dump(output, f, indent=2, ensure_ascii=False)
        dialog.destroy()


def main():
    win = App()
    win.show_all()
    Gtk.main()

if __name__ == "__main__":
    main()