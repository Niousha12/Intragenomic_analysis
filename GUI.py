import os
import pickle
import threading
import tkinter
import tkinter.messagebox
from functools import partial
from pathlib import Path
import glob
import random

import customtkinter
from tkinter import filedialog
from tkinter import ttk

# import mplcursors
import numpy as np
from PIL import Image
# from CGR import CGR, normalize, get_dist, read_fasta, reverse_complement

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib import pyplot as plt

from sklearn import cluster
from sklearn.decomposition import PCA

from chromosomes_holder import ChromosomesHolder
from constants import DISTANCE_METRICS_LIST, SCIENTIFIC_NAMES

# from dna_features_viewer import GraphicFeature, GraphicRecord

customtkinter.set_appearance_mode("System")  # Modes: "System" (standard), "Dark", "Light"
customtkinter.set_default_color_theme("blue")  # Themes: "blue" (standard), "green", "dark-blue"


class GUIDataStructure:
    def __init__(self):
        self.specie = customtkinter.StringVar()
        self.chromosome = customtkinter.StringVar()
        self.seq = ""
        self.start_seq = customtkinter.IntVar()
        self.start_txt = customtkinter.StringVar()
        self.end_seq = customtkinter.IntVar()
        self.end_txt = customtkinter.StringVar()
        self.annotation = customtkinter.StringVar()

    def invalidate_based_specie(self):
        self.chromosome.set("")
        self.start_seq.set(0)
        self.end_seq.set(0)
        self.annotation.set("")

    def invalidate_based_chromosome(self):
        self.start_seq.set(0)
        self.end_seq.set(0)
        self.annotation.set("")


class App(customtkinter.CTk):
    def __init__(self):
        super().__init__()

        # configure window TODO: find geometry window base on the screen size
        self.title("CGR GUI.py")
        self.geometry(f"{2300}x{1500}")
        self.header_font = ('Cambria', 14, 'bold')

        # configure grid layout (4x4)
        self.grid_columnconfigure(0, weight=1)
        self.grid_rowconfigure(0, weight=1)

        # create tabview
        self.tabview = customtkinter.CTkTabview(self)
        self.tabview.grid(padx=(20, 20), pady=(20, 20), sticky="nsew")
        tab_names = ["CGR Comparator", "Consecutive Windows", "Common Reference"]
        for tab_name in tab_names:
            self.tabview.add(tab_name)

        '''
            First Page (CGR Comparator tab)
            Configuring 
        '''
        self.tabview.tab(tab_names[0]).grid_columnconfigure(0, weight=1)
        self.tabview.tab(tab_names[0]).grid_columnconfigure(1, weight=10)
        self.tabview.tab(tab_names[0]).grid_rowconfigure(0, weight=1)
        self.tabview.tab(tab_names[0]).grid_rowconfigure(1, weight=4)

        # Frames
        self.t1_config_frame = customtkinter.CTkFrame(self.tabview.tab(tab_names[0]), corner_radius=20,
                                                      fg_color="#333333")
        self.t1_config_frame.grid(row=0, column=0, rowspan=2, sticky="nsew")
        self.t1_slider_frame = customtkinter.CTkFrame(self.tabview.tab(tab_names[0]), corner_radius=20,
                                                      fg_color="#333333")
        self.t1_slider_frame.grid(row=0, column=1, padx=(5, 5), pady=(5, 5), sticky="nsew")
        self.t1_display_frame = customtkinter.CTkFrame(self.tabview.tab(tab_names[0]), corner_radius=20,
                                                       fg_color="#707370")
        self.t1_display_frame.grid(row=1, column=1, padx=(5, 5), pady=(5, 5), sticky="nsew")

        # Designing the config frame (F1)
        for row in range(10):  # Increased row count to account for empty rows
            self.t1_config_frame.grid_rowconfigure(row, weight=1)

        # Creating frames for chromosome 1 and chromosome 2
        self.t1_chr_frame = customtkinter.CTkFrame(self.t1_config_frame, corner_radius=20)
        self.t1_chr_frame.grid(row=0, columnspan=2, padx=5, pady=5, sticky="ew")

        # Chromosomes Widget
        species_list = list(SCIENTIFIC_NAMES.keys())
        self.ds = {'1': GUIDataStructure(), '2': GUIDataStructure()}
        self.t1_species_combobox = {}
        self.t1_chr_combobox = {}
        for i in range(2):
            label = customtkinter.CTkLabel(self.t1_chr_frame, text=f"Chromosome {i + 1}: ", font=self.header_font)
            label.grid(row=i * 3, column=0, sticky="w", padx=10, pady=(5, 0))
            specie_label = customtkinter.CTkLabel(self.t1_chr_frame, text="Specie: ", font=self.header_font)
            specie_label.grid(row=(i * 3) + 1, column=0, sticky="w", padx=10)
            chr_label = customtkinter.CTkLabel(self.t1_chr_frame, text="Chromosome name: ", font=self.header_font)
            chr_label.grid(row=(i * 3) + 1, column=1, sticky="w", padx=10)
            self.t1_species_combobox[f"{i + 1}"] = customtkinter.CTkComboBox(self.t1_chr_frame, values=species_list,
                                                                             width=100,
                                                                             variable=self.ds[str(i + 1)].specie,
                                                                             command=partial(self.specie_change_event,
                                                                                             f"{i + 1}"))

            self.t1_species_combobox[f"{i + 1}"].grid(row=(i * 3) + 2, column=0, sticky="w", padx=10, pady=(0, 10))
            self.t1_species_combobox[f"{i + 1}"].set("")

            self.t1_chr_combobox[f"{i + 1}"] = customtkinter.CTkComboBox(self.t1_chr_frame, values=[],
                                                                         variable=self.ds[str(i + 1)].chromosome,
                                                                         width=100,
                                                                         command=partial(self.chromosome_change_event,
                                                                                         f"{i + 1}"))
            self.t1_chr_combobox[f"{i + 1}"].grid(row=(i * 3) + 2, column=1, sticky="w", padx=10, pady=(0, 10))
            self.t1_chr_combobox[f"{i + 1}"].set("")

        # k_mer combo box
        self.k_var = customtkinter.IntVar()
        k_mer_label = customtkinter.CTkLabel(self.t1_config_frame, text="k-mer: ", font=self.header_font)
        k_mer_label.grid(row=2, column=0, padx=10, pady=(10, 10))
        values_list = ["2", "3", "4", "5", "6", "7", "8", "9"]
        self.t1_k_mer_combobox = customtkinter.CTkComboBox(self.t1_config_frame, values=values_list, width=100,
                                                           state="normal", variable=self.k_var)
        self.t1_k_mer_combobox.grid(row=2, column=1, sticky="w", padx=10, pady=(10, 10))

        # Radio Button (Window size)
        # Frame for window size settings
        config_frame_color = self.t1_config_frame.cget("fg_color")
        self.t1_window_size_frame = customtkinter.CTkFrame(self.t1_config_frame, fg_color=config_frame_color)
        self.t1_window_size_frame.grid(row=3, columnspan=2, padx=10, pady=5, sticky="w")

        window_s_label = customtkinter.CTkLabel(self.t1_window_size_frame, text="Window Size:", font=self.header_font)
        window_s_label.grid(row=0, column=0, padx=10)
        self.window_s_toggle = tkinter.IntVar(value=0)
        self.t1_window_s_1 = customtkinter.CTkRadioButton(self.t1_window_size_frame, text="Variable",
                                                          variable=self.window_s_toggle,
                                                          value=0, command=self.window_size_toggle_event)
        self.t1_window_s_1.grid(row=1, column=0, padx=10, pady=5)
        self.t1_window_s_2 = customtkinter.CTkRadioButton(self.t1_window_size_frame, text="Fix",
                                                          variable=self.window_s_toggle,
                                                          value=1, command=self.window_size_toggle_event)
        self.t1_window_s_2.grid(row=1, column=1, padx=10, pady=5, sticky="w")
        self.window_s = tkinter.StringVar(value="")
        self.window_entry = customtkinter.CTkEntry(self.t1_window_size_frame, textvariable=self.window_s)
        self.window_entry.bind('<FocusOut>', partial(self.sequence_value_change, "0"))
        self.window_entry.bind('<Key-Return>', partial(self.sequence_value_change, "0"))
        self.window_entry.configure(state="disable")
        self.window_entry.grid(row=2, columnspan=2, padx=(10, 10), pady=(10, 10), sticky="ew")

        # reverse complement and random
        self.t1_seq1_rv = customtkinter.CTkFrame(self.t1_config_frame, fg_color=config_frame_color)
        self.t1_seq1_rv.grid(row=4, columnspan=2, padx=10, pady=5, sticky="ew")

        self.checkbox_RC = {}
        self.checkbox_Random = {}
        for i in range(2):
            seq_label = customtkinter.CTkLabel(self.t1_seq1_rv, text=f'Sequence {i + 1} :', font=self.header_font)
            seq_label.grid(row=(i * 2), column=0, padx=10, pady=(5, 0), sticky="w")

            self.checkbox_RC[str(i + 1)] = customtkinter.CTkCheckBox(master=self.t1_seq1_rv, text="Reverse Complement")
            self.checkbox_RC[str(i + 1)].grid(row=(i * 2), column=1, padx=10, pady=5, sticky="w")
            self.checkbox_Random[str(i + 1)] = customtkinter.CTkCheckBox(master=self.t1_seq1_rv, text="Shuffle")
            self.checkbox_Random[str(i + 1)].grid(row=(i * 2) + 1, column=1, padx=10, pady=5, sticky="w")

        # Distance metrics
        self.dist_metric = customtkinter.StringVar()
        dist_metric_label = customtkinter.CTkLabel(self.t1_config_frame, text="Distance Metric: ",
                                                   font=self.header_font)
        dist_metric_label.grid(row=6, column=0, padx=10, pady=(10, 10))
        self.t1_dist_metric_combobox = customtkinter.CTkComboBox(self.t1_config_frame, values=DISTANCE_METRICS_LIST,
                                                                 width=120, variable=self.dist_metric)
        self.t1_dist_metric_combobox.grid(row=6, column=1, sticky="w", padx=10, pady=(10, 10))
        self.t1_dist_metric_combobox.set("")

        # cgr/fcgr option
        self.fcgr = customtkinter.IntVar(value=1)
        switch = customtkinter.CTkSwitch(self.t1_config_frame, text=f"FCGR", variable=self.fcgr)
        switch.grid(row=7, columnspan=2, pady=(10, 10))

        # # plot button
        # cgr_img_base = customtkinter.CTkImage(Image.open(f"{self.path_dict['assets']}/background.png"),
        #                                       size=(1000, 360))
        # self.cgr_img_l = customtkinter.CTkLabel(self.t1_display_frame, image=cgr_img_base, text="")
        # self.cgr_img_l.grid(row=0, padx=20, sticky='nsew')
        # self.cgr_img = None
        # color_map = customtkinter.CTkImage(Image.open(f"{self.path_dict['assets']}/colormap.png"), size=(400, 50))
        # cgr_map = customtkinter.CTkLabel(self.t1_display_frame, image=color_map, text="")
        # cgr_map.grid(row=1, padx=20)
        t1_plot_button = customtkinter.CTkButton(self.t1_config_frame, text="Plot", width=120)  # command=self.plot
        t1_plot_button.grid(row=8, columnspan=2, pady=(10, 10))

        # Configuring First tab/second frame
        # Designing the slider frame (F2)
        for row in range(4):
            self.t1_config_frame.grid_rowconfigure(row, weight=1)
        for col in range(4):
            self.t1_config_frame.grid_columnconfigure(col, weight=1)

        # First Sequence scale
        _pad_size = (20, 0)
        self.t1_parts_name_combobox = {}
        self.start_seq_scale = {}
        self.end_seq_scale = {}
        self.t1_start_seq_entry = {}
        self.t1_end_seq_entry = {}
        for i in range(2):
            seq_label = customtkinter.CTkLabel(self.t1_slider_frame, text=f'Sequence {i + 1} :', font=self.header_font)
            seq_label.grid(row=(i * 2), column=0, padx=10, pady=_pad_size)

            # Sequence part names combo box
            self.t1_parts_name_combobox[str(i + 1)] = \
                customtkinter.CTkComboBox(self.t1_slider_frame, width=100, values=[],
                                          command=partial(self.annotation_change_event, str(i + 1)),
                                          variable=self.ds[str(i + 1)].annotation, state="disable")
            self.t1_parts_name_combobox[str(i + 1)].grid(row=(i * 2) + 1, column=0, padx=10, pady=_pad_size)

            seq_label_s = customtkinter.CTkLabel(self.t1_slider_frame, text='Start')
            seq_label_s.grid(row=(i * 2), column=1, padx=5, pady=_pad_size)
            seq_label_e = customtkinter.CTkLabel(self.t1_slider_frame, text='End')
            seq_label_e.grid(row=(i * 2) + 1, column=1, padx=5, pady=_pad_size)

            self.start_seq_scale[str(i + 1)] = ttk.Scale(self.t1_slider_frame, from_=0, to=len(self.ds[str(i + 1)].seq),
                                                         orient=customtkinter.HORIZONTAL, length=700,
                                                         variable=self.ds[str(i + 1)].start_seq,
                                                         command=partial(self.sequence_value_change, str(i + 1)))
            self.start_seq_scale[str(i + 1)].grid(row=(i * 2), column=2, pady=_pad_size)

            self.t1_start_seq_entry[str(i + 1)] = customtkinter.CTkEntry(self.t1_slider_frame,
                                                                         textvariable=self.ds[str(i + 1)].start_txt)
            self.t1_start_seq_entry[str(i + 1)].bind('<FocusOut>', partial(self.sequence_value_change, "3"))
            self.t1_start_seq_entry[str(i + 1)].bind('<Key-Return>', partial(self.sequence_value_change, "3"))
            self.t1_start_seq_entry[str(i + 1)].grid(row=(i * 2), column=3, padx=5, pady=_pad_size)
            seq_s_e_label = customtkinter.CTkLabel(self.t1_slider_frame, text='bp')
            seq_s_e_label.grid(row=(i * 2), column=4, pady=_pad_size)

            self.end_seq_scale[str(i + 1)] = ttk.Scale(self.t1_slider_frame, from_=0,
                                                       to=self.ds[str(i + 1)].end_seq.get(),
                                                       orient=customtkinter.HORIZONTAL, length=700,
                                                       variable=self.ds[str(i + 1)].end_seq,
                                                       command=partial(self.sequence_value_change, str(i + 1)))
            self.end_seq_scale[str(i + 1)].grid(row=(i * 2) + 1, column=2, pady=_pad_size)
            self.t1_end_seq_entry[str(i + 1)] = customtkinter.CTkEntry(self.t1_slider_frame,
                                                                       textvariable=self.ds[str(i + 1)].end_txt)
            self.t1_end_seq_entry[str(i + 1)].bind('<FocusOut>', partial(self.sequence_value_change, "3"))
            self.t1_end_seq_entry[str(i + 1)].bind('<Key-Return>', partial(self.sequence_value_change, "3"))
            self.t1_end_seq_entry[str(i + 1)].grid(row=(i * 2) + 1, column=3, padx=5, pady=_pad_size)
            seq_s_e_label_d = customtkinter.CTkLabel(self.t1_slider_frame, text='bp')
            seq_s_e_label_d.grid(row=(i * 2) + 1, column=4, pady=_pad_size)

    def sync_text_vars(self, sender, keep_annotation=False):
        # for key, value in self.ds.items():
        self.ds[sender].start_txt.set(f"{self.ds[sender].start_seq.get()}")
        self.ds[sender].end_txt.set(f"{self.ds[sender].end_seq.get()}")
        if keep_annotation is False:
            self.ds[sender].annotation.set("")

    def reverse_sync_text_vars(self, sender):
        # for key, value in self.ds.items():
        self.ds[sender].start_seq.set(int(self.ds[sender].start_txt.get()))
        self.ds[sender].end_seq.set(int(self.ds[sender].end_txt.get()))
        self.ds[sender].annotation.set("")

    # '''events'''
    def specie_change_event(self, sender, value):
        specie = self.ds[sender].specie.get()
        self.t1_chr_combobox[sender].configure(values=ChromosomesHolder(specie).get_all_chromosomes_name())
        self.ds[sender].invalidate_based_specie()

    def chromosome_change_event(self, sender, value):
        specie = self.ds[sender].specie.get()
        chromosome = self.ds[sender].chromosome.get()
        # set its sequence
        self.ds[sender].seq = ChromosomesHolder(specie).get_chromosome_sequence(chromosome)
        # set the annotation combobox
        self.t1_parts_name_combobox[sender].configure(state="normal")
        self.t1_parts_name_combobox[sender].configure(values=ChromosomesHolder(specie).cytobands[chromosome])
        self.ds[sender].invalidate_based_chromosome()
        # set start and end
        self.ds[sender].end_seq.set(len(self.ds[sender].seq))
        if len(self.ds[sender].seq) > 0:
            self.start_seq_scale[sender].configure(state="normal")
            self.end_seq_scale[sender].configure(state="normal")
        else:
            self.start_seq_scale[sender].configure(state="disable")
            self.end_seq_scale[sender].configure(state="disable")
        self.start_seq_scale[sender].configure(to=len(self.ds[sender].seq))
        self.end_seq_scale[sender].configure(to=len(self.ds[sender].seq))

        self.sync_text_vars(sender)

    def annotation_change_event(self, sender, value):
        annotation = self.ds[sender].annotation.get()
        chromosome_name = self.ds[sender].chromosome.get()
        annotation_info = ChromosomesHolder(self.ds[sender].specie.get()).cytobands[chromosome_name][annotation]
        self.ds[sender].start_seq.set(annotation_info.start)
        self.ds[sender].end_seq.set(annotation_info.end)
        self.sync_text_vars(sender, keep_annotation=True)

    def window_size_toggle_event(self):
        if self.window_s_toggle.get() == 0:
            self.window_s.set("")
            self.window_entry.configure(state="disable")

            # end scales
            for key, value in self.end_seq_scale.items():
                value.configure(state="normal")
            for key, value in self.t1_end_seq_entry.items():
                value.configure(state="normal")
        else:
            self.window_entry.configure(state="normal")
            self.window_s.set("500000")

            # end scales
            for key, value in self.end_seq_scale.items():
                value.configure(state="disable")
            for key, value in self.t1_end_seq_entry.items():
                value.configure(state="disable")
            for key, value in self.ds.items():
                self.ds[key].end_seq.set(self.ds[key].start_seq.get() + int(self.window_s.get()))

        for key, value in self.ds.items():
            self.sync_text_vars(key)

    def sequence_value_change(self, sender, value):
        if sender == "0":  # Window size changed
            for key, value in self.ds.items():
                self.ds[key].end_seq.set(self.ds[key].start_seq.get() + int(self.window_s.get()))
        elif sender in ["1", "2"]:  # Scale changed
            if self.window_s_toggle.get() == 1:
                self.ds[sender].end_seq.set(self.ds[sender].start_seq.get() + int(self.window_s.get()))
        elif sender in ["3"]:
            for key, value in self.ds.items():
                self.reverse_sync_text_vars(key)
            if self.window_s_toggle.get() == 1:
                for key, value in self.ds.items():
                    self.ds[key].end_seq.set(self.ds[key].start_seq.get() + int(self.window_s.get()))

        for key, value in self.ds.items():
            self.sync_text_vars(key)

        # if self.window_s_toggle.get() == 1:
        #     self.end_var.set(self.begin_var.get() + int(self.window_s.get()) / self._scale)
        #     self.end_var_2.set(self.begin_var_2.get() + int(self.window_s.get()) / self._scale)

    # def entry_text_changed(self):
    #     print("I am called")
    #     # self.ds[sender].start_seq.set(int(self.t1_start_seq_entry[sender].get()))

    # class ChrAnnotationRecord:
    #     def __init__(self, chr_name, start, end, name, display_name, group):
    #         self.name = name
    #         self.chr_name = chr_name
    #         self.start = start
    #         self.end = end
    #         self.display_name = display_name
    #         self.group = group
    #
    #     def __repr__(self):
    #         return f"{self.display_name}, {self.start}, {self.end}"

    # def _init_knowledge_dicts(self):
    #     # reading telomere
    #     file_path = "assets/chm13v2.0_telomere.bed"
    #     switch = 0
    #     with open(file_path) as file:
    #         for line in file:
    #             line = line.strip()
    #             parts = line.split("\t")
    #
    #             self.knowledge[f"tel_{switch + 1}_{parts[0]}"] = \
    #                 self.ChrAnnotationRecord(chr_name=parts[0], start=int(parts[1]), end=int(parts[2]),
    #                                          name=f"tel_{switch + 1}_{parts[0]}", display_name=f"Telomere {switch + 1}",
    #                                          group="telomeres")
    #
    #             switch = 1 - switch
    #
    #     # reading centromere
    #     file_path = "assets/chm13v2.0_censat_v2.0.bed"
    #     chr_name = None
    #     last_chr = None
    #     last_start = None
    #     last_end = None
    #     with open(file_path) as file:
    #         for index, line in enumerate(file):
    #             if index == 0:
    #                 continue
    #             line = line.strip()
    #             parts = line.split("\t")
    #
    #             if chr_name != parts[0]:
    #                 if not last_end is None:
    #                     self.knowledge[f"cent_{chr_name}"] = \
    #                         self.ChrAnnotationRecord(chr_name=chr_name, start=last_start, end=last_end,
    #                                                  name=f"cent_{chr_name}",
    #                                                  display_name=f"Centromere",
    #                                                  group="centromeres")
    #
    #                 last_start = int(parts[1])
    #                 chr_name = parts[0]
    #
    #             last_end = int(parts[2])
    #
    #     # Cytobands
    #     file_path = "assets/chm13v2.0_cytobands_allchrs.bed"
    #     with open(file_path) as file:
    #         for line in file:
    #             line = line.strip()
    #             parts = line.split("\t")
    #             self.knowledge[f"cytoband_{parts[0]}_{parts[3]}"] = \
    #                 self.ChrAnnotationRecord(chr_name=parts[0], start=int(parts[1]), end=int(parts[2]),
    #                                          name=f"cytoband_{parts[0]}_{parts[3]}",
    #                                          display_name=f"{parts[3]}",
    #                                          group="cytobands")
    #
    #     _, coi = read_fasta("assets/barcode.fna")
    #     self.barcode = coi
    #
    # def clear_temp(self):
    #     root = Path(f"{self.path_dict['temp']}")
    #     if root.exists():
    #         # file_names = [f for f in listdir(self.path_dict["temp"]) if isfile(join(self.path_dict["temp"], f))]
    #         file_names = glob.glob(f"{self.path_dict['temp']}" + '/**/*.*', recursive=True)
    #         for file in file_names:
    #             path = Path(file)
    #             path.unlink()
    #
    # def browse_file(self):
    #     file_path = filedialog.askopenfilename()
    #     self.sequence = ""
    #     with open(file_path) as file:
    #         for line in file:
    #             line = line.strip()
    #             if line.startswith(">"):
    #                 self.file_name = line.split(",")[0].split(" ")[-1]
    #             else:
    #                 self.sequence += line
    #
    #     # green color
    #     self.status_text_box.configure(text=f"Chr {self.file_name} uploaded!", text_color="green")
    #     self.t2_chr_text_box.configure(text=f"Chromosome {self.file_name}", text_color="green")
    #     self.t3_chr_text_box.configure(text=f"Chromosome {self.file_name}", text_color="green")
    #     self.t1_k_mer_combobox.configure(state="normal")
    #     self.t2_k_mer_combobox.configure(state="normal")
    #     self.t3_k_mer_combobox.configure(state="normal")
    #     # self.t1_dist_metric_combobox.configure(state="normal")
    #
    #     self.end_scale = len(self.sequence) / self._scale
    #     self.first_seq_s_scale.configure(to=self.end_scale)
    #     self.first_seq_e_scale.configure(to=self.end_scale)
    #     self.sec_seq_s_scale.configure(to=self.end_scale)
    #     self.sec_seq_e_scale.configure(to=self.end_scale)
    #     self.t3_ref_end.set(len(self.sequence) / self._scale)
    #     self.begin_var.set(0.0)
    #     self.end_var.set(0.0)
    #     self.begin_var_2.set(0.0)
    #     self.end_var_2.set(0.0)
    #
    #     # filling the combo box
    #     knowledge_keys = list(self.knowledge.keys())
    #     self.different_parts_list = []
    #     self.different_parts_list.append("mt barcode")
    #     for key in knowledge_keys:
    #         if key == 'mt_barcode':
    #             pass
    #         else:
    #             annotation_info = self.knowledge[key]
    #             if self.file_name == annotation_info.chr_name.split('r')[1]:
    #                 self.different_parts_list.append(annotation_info.display_name)
    #
    #     self.t1_seq_1_parts_name_combobox.configure(state="normal")
    #     self.t1_seq_1_parts_name_combobox.set("")
    #     self.t1_seq_2_parts_name_combobox.configure(state="normal")
    #     self.t1_seq_2_parts_name_combobox.set("")
    #     self.t1_seq_1_parts_name_combobox.configure(values=self.different_parts_list)
    #     self.t1_seq_2_parts_name_combobox.configure(values=self.different_parts_list)
    #     self.t3_annotation_combobox.configure(state="normal")
    #     self.t3_annotation_combobox.configure(values=self.different_parts_list)
    #
    #     self.cgr_images_clustering = []
    #     self.band_info_clustering = []
    #
    # def window_size_toggle_event(self):
    #     if self.window_s_toggle.get() == 0:
    #         self.window_s.set("")
    #         self.window_entry.configure(state="disable")
    #
    #         # end scales
    #         self.first_seq_e_scale.configure(state="normal")
    #         self.first_seq_e.configure(state="normal")
    #         # self.begin_var.set(0)
    #         # self.end_var.set(0)
    #         self.sec_seq_e_scale.configure(state="normal")
    #         self.sec_seq_e.configure(state="normal")
    #         # self.begin_var_2.set(0)
    #         # self.end_var_2.set(0)
    #     else:
    #         self.window_entry.configure(state="normal")
    #         self.window_s.set("2000000")
    #
    #         # end scales
    #         self.first_seq_e_scale.configure(state="disable")
    #         self.end_var.set(self.begin_var.get() + int(self.window_s.get()) / self._scale)
    #         self.first_seq_e.configure(state="disable")
    #         self.sec_seq_e_scale.configure(state="disable")
    #         self.end_var_2.set(self.begin_var_2.get() + int(self.window_s.get()) / self._scale)
    #         self.sec_seq_e.configure(state="disable")
    #
    # def annotation_change_event_seq1(self, value):
    #     if self.t1_seq_1_annotation_value.get() == "mt barcode":
    #         self.begin_var.set(0 / self._scale)
    #         self.end_var.set(0 / self._scale)
    #     else:
    #         for key, value in self.knowledge.items():
    #             if value.display_name == self.t1_seq_1_annotation_value.get() and self.file_name == \
    #                     value.chr_name.split("r")[1]:
    #                 self.begin_var.set(value.start / self._scale)
    #                 self.end_var.set(value.end / self._scale)
    #                 break
    #
    # def annotation_change_event_seq2(self, value):
    #     if self.t1_seq_2_annotation_value.get() == "mt barcode":
    #         self.begin_var_2.set(0 / self._scale)
    #         self.end_var_2.set(0 / self._scale)
    #     else:
    #         for key, value in self.knowledge.items():
    #             if value.display_name == self.t1_seq_2_annotation_value.get() and self.file_name == \
    #                     value.chr_name.split("r")[1]:
    #                 self.begin_var_2.set(value.start / self._scale)
    #                 self.end_var_2.set(value.end / self._scale)
    #                 break
    #
    # def annotation_change_event_ref(self, value):
    #     if self.t3_annotation_value.get() == "mt barcode":
    #         self.t3_ref_begin.set(0 / self._scale)
    #         self.t3_ref_end.set(0 / self._scale)
    #     else:
    #         for key, value in self.knowledge.items():
    #             if value.display_name == self.t3_annotation_value.get() and \
    #                     self.file_name == value.chr_name.split("r")[1]:
    #                 self.t3_ref_begin.set(value.start / self._scale)
    #                 self.t3_ref_end.set(value.end / self._scale)
    #                 break

    # def move_previous(self, value):
    #     if self.t2_pic_num.get() > 0:
    #         self.t2_pic_num.set(self.t2_pic_num.get() - 1)
    #         self._change_image(None)
    #
    # def t3_move_previous(self, value):
    #     if self.t3_pic_num.get() > 0:
    #         self.t3_pic_num.set(self.t3_pic_num.get() - 1)
    #         self.t3_change_image(None)
    #
    # def move_next(self, value):
    #     if self.t2_pic_num.get() < len(self.cgr_distance_history) - 1:
    #         self.t2_pic_num.set(self.t2_pic_num.get() + 1)
    #         self._change_image(None)
    #
    # def t3_move_next(self, value):
    #     if self.t3_pic_num.get() < len(self.t3_cgr_distance_history) - 1:
    #         self.t3_pic_num.set(self.t3_pic_num.get() + 1)
    #         self.t3_change_image(None)
    #
    # def plot(self):
    #     begin_seq = int(self.begin_var.get() * self._scale)
    #     end_seq_1 = int(self.end_var.get() * self._scale)
    #     sec_seq = int(self.begin_var_2.get() * self._scale)
    #     end_seq_2 = int(self.end_var_2.get() * self._scale)
    #     kmer = int(self.k_var.get())
    #
    #     if self.t1_seq_1_annotation_value.get() == "mt barcode":
    #         modified = self.barcode
    #     else:
    #         modified = self.sequence[begin_seq:end_seq_1]
    #     if self.checkbox_RC_seq1.get():
    #         modified = reverse_complement(modified)
    #     if self.checkbox_Random_seq1.get():
    #         modified = list(modified)
    #         random.shuffle(modified)
    #         modified = ''.join(modified)
    #     cgr1 = CGR(modified, kmer)
    #
    #     if self.t1_seq_2_annotation_value.get() == "mt barcode":
    #         modified_2 = self.barcode
    #     else:
    #         modified_2 = self.sequence[sec_seq:end_seq_2]
    #     if self.checkbox_RC_seq2.get():
    #         modified_2 = reverse_complement(modified_2)
    #     if self.checkbox_Random_seq2.get():
    #         modified_2 = list(modified_2)
    #         random.shuffle(modified_2)
    #         modified_2 = ''.join(modified_2)
    #     cgr2 = CGR(modified_2, kmer)
    #
    #     if self.fcgr.get() == 1:
    #         img_1 = cgr1.frequency_chaos_game_representation()
    #         img_2 = cgr2.frequency_chaos_game_representation()
    #         img_1 = normalize(img_1)
    #         img_2 = normalize(img_2)
    #     else:
    #         img_1 = cgr1.chaos_game_representation()
    #         img_2 = cgr2.chaos_game_representation()
    #
    #     diff = img_2 - img_1
    #
    #     dist = get_dist(img_1, img_2, dist_m=self.dist_metric.get())
    #
    #     fig, (ax1, ax2, ax3) = plt.subplots(1, 3)
    #     extent = 0, 1, 0, 1
    #     norm = plt.Normalize(-0.5, 0.5)
    #
    #     # plot the data on the subplots
    #     ax1.imshow(img_1, cmap='Blues', extent=extent)
    #     ax1.tick_params(right=False, labelbottom=False, bottom=False)
    #     ax1.set_title(f'{round(begin_seq / self._scale, 2)}M_{round(end_seq_1 / self._scale, 2)}M')
    #
    #     ax2.imshow(diff, cmap='RdBu', norm=norm, extent=extent)
    #     ax2.tick_params(left=False, right=False, labelleft=False, labelbottom=False, bottom=False)
    #     ax2.set_title(f'distance = {round(dist, 4)}')
    #
    #     ax3.imshow(img_2, cmap='Blues', extent=extent)
    #     ax3.tick_params(left=False, right=False, labelleft=False, labelbottom=False, bottom=False)
    #     ax3.set_title(f'{round(sec_seq / self._scale, 2)}M_{round(end_seq_2 / self._scale, 2)}M')
    #
    #     plt.savefig(f"{self.path_dict['temp']}/temp_1.png", bbox_inches='tight', transparent=True)
    #
    #     self.cgr_img = customtkinter.CTkImage(Image.open(f"{self.path_dict['temp']}/temp_1.png"), size=(1000, 360))
    #     self.cgr_img_l.configure(image=self.cgr_img)
    #
    # def run_consecutive(self, event):
    #     global foo_thread
    #     foo_thread = threading.Thread(target=self.t2_run)
    #     foo_thread.daemon = True
    #     foo_thread.start()
    #     self.after(20, self.check_thread)
    #
    # def _generate_pickled_image(self, image_index, mode):
    #     if mode == "consecutive":
    #         with open(f"{self.path_dict['temp']}/consecutive/pickle/{image_index}.pkl", 'rb') as handle:
    #             dictionary = pickle.load(handle)
    #
    #         cgr1 = dictionary["cgr1"]
    #         im1 = dictionary["im1"]
    #         cgr2 = dictionary["cgr2"]
    #         im2 = dictionary["im2"]
    #         diff = dictionary["diff"]
    #         dist = dictionary["dist"]
    #         b1 = dictionary["b1"]
    #         e1 = dictionary["e1"]
    #         b2 = dictionary["b2"]
    #         e2 = dictionary["e2"]
    #
    #         fig, (ax1, ax2, ax3) = plt.subplots(1, 3)
    #         extent = 0, 1, 0, 1
    #         norm = plt.Normalize(-0.5, 0.5)
    #
    #         ax1.imshow(im1, cmap='Blues', extent=extent)
    #         ax1.tick_params(right=False, labelbottom=False, bottom=False)
    #         ax1.set_title(f'{round(b1 / self._scale, 2)}M_{round(e1 / self._scale, 2)}M')
    #
    #         ax2.imshow(diff, cmap='RdBu', norm=norm, extent=extent)
    #         ax2.tick_params(left=False, right=False, labelleft=False, labelbottom=False, bottom=False)
    #         ax2.set_title(f'distance = {round(dist, 4)}')
    #
    #         ax3.imshow(im2, cmap='Blues', extent=extent)
    #         ax3.tick_params(left=False, right=False, labelleft=False, labelbottom=False, bottom=False)
    #         ax3.set_title(f'{round(b2 / self._scale, 2)}M_{round(e2 / self._scale, 2)}M')
    #
    #         path = f"{self.path_dict['temp']}/consecutive/{self.file_name}"
    #         if not os.path.exists(path):
    #             os.makedirs(path)
    #         plt.savefig(f"{path}/{image_index}.png", bbox_inches='tight', transparent=True)
    #     elif mode == "common_ref":
    #         with open(f"{self.path_dict['temp']}/common_ref/pickle/{image_index}.pkl", 'rb') as handle:
    #             dictionary = pickle.load(handle)
    #
    #         cgr2 = dictionary["cgr2"]
    #         im2 = dictionary["im2"]
    #         diff = dictionary["diff"]
    #         dist = dictionary["dist"]
    #         b1 = dictionary["b1"]
    #         e1 = dictionary["e1"]
    #
    #         fig, (ax1, ax2) = plt.subplots(1, 2)
    #         extent = 0, 1, 0, 1
    #         norm = plt.Normalize(-0.5, 0.5)
    #
    #         ax1.imshow(diff, cmap='RdBu', norm=norm, extent=extent)
    #         ax1.tick_params(left=False, right=False, labelleft=False, labelbottom=False, bottom=False)
    #         ax1.set_title(f'distance = {round(dist, 4)}')
    #
    #         ax2.imshow(im2, cmap='Blues', extent=extent)
    #         ax2.tick_params(left=False, right=False, labelleft=False, labelbottom=False, bottom=False)
    #         ax2.set_title(f'{round(b1 / self._scale, 2)}M_{round(e1 / self._scale, 2)}M')
    #
    #         path = f"{self.path_dict['temp']}/common_ref/{self.file_name}"
    #         if not os.path.exists(path):
    #             os.makedirs(path)
    #         plt.savefig(f"{path}/{image_index}.png", bbox_inches='tight', transparent=True)
    #
    # def _change_image(self, value):
    #     self._generate_pickled_image(self.t2_pic_num.get(), mode="consecutive")
    #     t2_cgr_img_r = customtkinter.CTkImage(
    #         Image.open(f"{self.path_dict['temp']}/consecutive/{self.file_name}/{self.t2_pic_num.get()}.png"),
    #         size=(800, 300))
    #     self.t2_cgr_img_l.configure(image=t2_cgr_img_r)
    #
    #     # plot bar again and change color of a bar
    #     fig, ax1 = plt.subplots(figsize=(10.7, 3))
    #     x = np.arange(len(self.cgr_distance_history))
    #     y = np.asarray(self.cgr_distance_history)
    #     mask1 = x == self.t2_pic_num.get()
    #     mask2 = x != self.t2_pic_num.get()
    #     ax1.bar(x[mask1], y[mask1], color='red')
    #     ax1.bar(x[mask2], y[mask2], color='blue')
    #
    #     canvas = FigureCanvasTkAgg(fig, self.t2_plot_frame)
    #     canvas.draw()
    #     canvas.get_tk_widget().grid(row=0, padx=15, pady=5)
    #
    # def t3_change_image(self, event):
    #     self._generate_pickled_image(self.t3_pic_num.get(), mode="common_ref")
    #     t3_cgr_img_t = customtkinter.CTkImage(
    #         Image.open(f"{self.path_dict['temp']}/common_ref/{self.file_name}/{self.t3_pic_num.get()}.png"),
    #         size=(600, 300))
    #     self.t3_cgr_img_l.configure(image=t3_cgr_img_t)
    #
    #     # plot bar again and change color of a bar
    #     fig, ax1 = plt.subplots(figsize=(10.7, 3))
    #     x = np.arange(len(self.t3_cgr_distance_history))
    #     y = np.asarray(self.t3_cgr_distance_history)
    #     mask1 = x == self.t3_pic_num.get()
    #     mask2 = x != self.t3_pic_num.get()
    #     ax1.bar(x[mask1], y[mask1], color='red')
    #     ax1.bar(x[mask2], y[mask2], color='blue')
    #
    #     canvas = FigureCanvasTkAgg(fig, self.t3_plot_frame)
    #     canvas.draw()
    #     canvas.get_tk_widget().grid(row=0, padx=15, pady=5)
    #
    # def check_thread(self):
    #     if foo_thread.is_alive():
    #         self.after(20, self.check_thread)
    #     else:
    #         self.t2_scale.configure(to=int(len(self.cgr_distance_history) - 1))
    #
    #         self._generate_pickled_image(0, mode="consecutive")
    #
    #         t2_cgr_img = customtkinter.CTkImage(
    #             Image.open(f"{self.path_dict['temp']}/consecutive/{self.file_name}/0.png"), size=(800, 300))
    #         self.t2_cgr_img_l.configure(image=t2_cgr_img)
    #
    #         # Display bar plot for cgr_distance_hist
    #         fig, ax1 = plt.subplots(figsize=(10.7, 3))
    #         x = np.arange(len(self.cgr_distance_history))
    #         y = np.asarray(self.cgr_distance_history)
    #         mask1 = x == 0
    #         mask2 = x != 0
    #         ax1.bar(x[mask1], y[mask1], color='red')
    #         ax1.bar(x[mask2], y[mask2], color='blue')
    #
    #         plt.savefig(
    #             f"{self.path_dict['outputs']}/consecutive_{self.file_name}_{self.window_s.get()}_"
    #             f"{self.dist_metric.get()}.png", bbox_inches='tight')
    #
    #         canvas = FigureCanvasTkAgg(fig, self.t2_plot_frame)
    #         canvas.draw()
    #         canvas.get_tk_widget().grid(row=0, padx=15, pady=5)
    #
    # def t2_run(self):
    #     self.cgr_distance_history = []
    #     self.step_length = np.floor(len(self.sequence) / int(self.window_s.get()))
    #     self.progress_bar.set(0)
    #     self.t2_pic_num.set(0)
    #     for i in range(int(self.step_length) - 1):
    #         dictionary = {}
    #         self.progress_bar.set((i + 2) / int(self.step_length))
    #         b1 = i * int(self.window_s.get())
    #         e1 = (i + 1) * int(self.window_s.get())
    #         b2 = (i + 1) * int(self.window_s.get())
    #         e2 = (i + 2) * int(self.window_s.get())
    #
    #         cgr1 = CGR(self.sequence[b1:e1], self.k_var.get())
    #         cgr2 = CGR(self.sequence[b2:e2], self.k_var.get())
    #
    #         if self.fcgr.get() == 1:
    #             im1 = cgr1.frequency_chaos_game_representation()
    #             im2 = cgr2.frequency_chaos_game_representation()
    #             im1 = normalize(im1)
    #             im2 = normalize(im2)
    #         else:
    #             im1 = cgr1.chaos_game_representation()
    #             im2 = cgr2.chaos_game_representation()
    #
    #         # im1 = cgr1.frequency_chaos_game_representation()
    #         # im2 = cgr2.frequency_chaos_game_representation()
    #
    #         diff = im2 - im1
    #
    #         dist = get_dist(im1, im2, dist_m=self.dist_metric.get())
    #
    #         self.cgr_distance_history.append(dist)
    #
    #         dictionary["cgr1"] = cgr1
    #         dictionary["im1"] = im1
    #         dictionary["cgr2"] = cgr2
    #         dictionary["im2"] = im2
    #         dictionary["diff"] = diff
    #         dictionary["dist"] = dist
    #         dictionary["b1"] = b1
    #         dictionary["e1"] = e1
    #         dictionary["b2"] = b2
    #         dictionary["e2"] = e2
    #
    #         path = f"{self.path_dict['temp']}/consecutive/pickle"
    #         if not os.path.exists(path):
    #             os.makedirs(path)
    #         with open(f"{path}/{i}.pkl", 'wb') as f:
    #             pickle.dump(dictionary, f)
    #
    # def run_common_ref(self, event):
    #     global foo_thread_2
    #     foo_thread_2 = threading.Thread(target=self.t3_run)
    #     foo_thread_2.daemon = True
    #     foo_thread_2.start()
    #     self.after(20, self.check_thread_2)
    #
    # def check_thread_2(self):
    #     if foo_thread_2.is_alive():
    #         self.after(20, self.check_thread_2)
    #     else:
    #         self.t3_scale.configure(to=int(len(self.t3_cgr_distance_history) - 1))
    #
    #         # generate the ref picture and display
    #         with open(f"{self.path_dict['temp']}/common_ref/pickle/ref.pkl", 'rb') as handle:
    #             ref_dict = pickle.load(handle)
    #
    #         cgr1 = ref_dict["cgr1"]
    #         im1 = ref_dict["im1"]
    #
    #         fig, (ax1) = plt.subplots(1, 1)
    #         extent = 0, 1, 0, 1
    #         norm = plt.Normalize(-0.5, 0.5)
    #
    #         ax1.imshow(im1, cmap='Blues', extent=extent)
    #         ax1.tick_params(left=False, right=False, labelleft=False, labelbottom=False, bottom=False)
    #         ax1.set_title(f'Reference')
    #
    #         path = f"{self.path_dict['temp']}/common_ref/{self.file_name}"
    #         if not os.path.exists(path):
    #             os.makedirs(path)
    #         plt.savefig(f"{path}/ref.png", bbox_inches='tight', transparent=True)
    #
    #         t3_ref_img = customtkinter.CTkImage(
    #             Image.open(f"{self.path_dict['temp']}/common_ref/{self.file_name}/ref.png"), size=(300, 300))
    #         self.t3_ref_img_l.configure(image=t3_ref_img)
    #
    #         # generate other pictures
    #         self._generate_pickled_image(0, mode="common_ref")
    #
    #         t3_cgr_img = customtkinter.CTkImage(
    #             Image.open(f"{self.path_dict['temp']}/common_ref/{self.file_name}/0.png"), size=(600, 300))
    #         self.t3_cgr_img_l.configure(image=t3_cgr_img)
    #
    #         # Display bar plot for cgr_distance_hist
    #         fig, ax1 = plt.subplots(figsize=(10.7, 3))
    #         x = np.arange(len(self.t3_cgr_distance_history))
    #         y = np.asarray(self.t3_cgr_distance_history)
    #         mask1 = x == 0
    #         mask2 = x != 0
    #         ax1.bar(x[mask1], y[mask1], color='red')
    #         ax1.bar(x[mask2], y[mask2], color='blue')
    #
    #         plt.savefig(
    #             f"{self.path_dict['outputs']}/reference_{self.file_name}_{self.t3_sliding_size.get()}_{self.dist_metric.get()}.png",
    #             bbox_inches='tight')
    #
    #         canvas = FigureCanvasTkAgg(fig, self.t3_plot_frame)
    #         canvas.draw()
    #         canvas.get_tk_widget().grid(row=0, padx=15, pady=5)
    #
    # def t3_run(self):
    #     self.t3_cgr_distance_history = []
    #     self.t3_step_length = np.floor(len(self.sequence) / int(self.t3_sliding_size.get()))
    #     self.t3_progress_bar.set(0)
    #     self.t3_pic_num.set(0)
    #
    #     if self.t3_annotation_value.get() == "mt barcode":
    #         cgr1 = CGR(self.barcode, self.k_var.get())
    #     else:
    #         b = int(self.t3_ref_begin.get() * self._scale)
    #         e = int(self.t3_ref_end.get() * self._scale)
    #
    #         cgr1 = CGR(self.sequence[b:e], self.k_var.get())
    #     self.t3_progress_bar.set(1 / (int(self.t3_step_length) + 2))
    #     if self.fcgr.get() == 1:
    #         cgr1_im = cgr1.frequency_chaos_game_representation()
    #         cgr1_im = normalize(cgr1_im)
    #     else:
    #         cgr1_im = cgr1.chaos_game_representation()
    #
    #     self.t3_progress_bar.set(2 / (int(self.t3_step_length) + 2))
    #
    #     ref_dict = {"cgr1": cgr1, "im1": cgr1_im}
    #
    #     path = f"{self.path_dict['temp']}/common_ref/pickle"
    #     if not os.path.exists(path):
    #         os.makedirs(path)
    #     with open(f"{path}/ref.pkl", 'wb') as f:
    #         pickle.dump(ref_dict, f)
    #
    #     for i in range(int(self.t3_step_length)):
    #         dictionary = {}
    #         self.t3_progress_bar.set((i + 3) / (int(self.t3_step_length) + 2))
    #         b1 = i * int(self.t3_sliding_size.get())
    #         e1 = (i + 1) * int(self.t3_sliding_size.get())
    #
    #         cgr2 = CGR(self.sequence[b1:e1], self.k_var.get())
    #         if self.fcgr.get() == 1:
    #             cgr2_im = cgr2.frequency_chaos_game_representation()
    #             cgr2_im = normalize(cgr2_im)
    #         else:
    #             cgr2_im = cgr2.chaos_game_representation()
    #
    #         diff = cgr2_im - cgr1_im
    #
    #         dist = get_dist(cgr1_im, cgr2_im, dist_m=self.dist_metric.get())
    #
    #         self.t3_cgr_distance_history.append(dist)
    #
    #         dictionary["cgr2"] = cgr2
    #         dictionary["im2"] = cgr2_im
    #         dictionary["diff"] = diff
    #         dictionary["dist"] = dist
    #         dictionary["b1"] = b1
    #         dictionary["e1"] = e1
    #
    #         with open(f"{path}/{i}.pkl", 'wb') as f:
    #             pickle.dump(dictionary, f)
    #
    # def run_clustering(self, event):
    #     global foo_thread_3
    #     foo_thread_3 = threading.Thread(target=self.run_chr_clustering)
    #     foo_thread_3.daemon = True
    #     foo_thread_3.start()
    #     self.after(20, self.check_thread_3)
    #
    # def check_thread_3(self):
    #     if foo_thread_3.is_alive():
    #         self.after(20, self.check_thread_3)
    #     else:
    #         cgr_images = np.asarray(self.cgr_images_clustering)
    #         cgr_images = cgr_images.reshape(cgr_images.shape[0], -1)
    #
    #         km = cluster.KMeans(n_clusters=3).fit(cgr_images)
    #
    #         pca = PCA(n_components=2)
    #         principal_components = pca.fit_transform(cgr_images)
    #
    #         fig, ax = plt.subplots(figsize=(3, 3))
    #         colors_dict = {0: "#ffd700", 1: "#ffcccc", 2: "#cffccc"}
    #         color_list = [colors_dict[l] for l in km.labels_]
    #         scatter = ax.scatter(principal_components[:, 0], principal_components[:, 1], c=color_list)
    #
    #         # Add labels to each point
    #         cursor = mplcursors.cursor(scatter)
    #         labels = []
    #         for band in self.band_info_clustering:
    #             labels.append(band.display_name)
    #
    #         @cursor.connect("add")
    #         def on_add(sel):
    #             sel.annotation.set_text(labels[sel.target.index])
    #
    #         canvas = FigureCanvasTkAgg(fig, self.t4_plot_frame)
    #         canvas.draw()
    #         canvas.get_tk_widget().grid(row=0, padx=(50, 0), pady=5)
    #
    #         # gene clustering
    #         features = []
    #         for i, info in enumerate(self.band_info_clustering):
    #             gene_graphic = GraphicFeature(start=int(info.start) / self._scale, end=int(info.end) / self._scale,
    #                                           color=colors_dict[km.labels_[i]], label=info.display_name)
    #             features.append(gene_graphic)
    #
    #         record = GraphicRecord(sequence_length=len(self.sequence) / self._scale, features=features)
    #         # plt.clf()
    #         record.plot()
    #         plt.savefig(f"{self.path_dict['temp']}/gene_cluster.png", bbox_inches='tight', transparent=True)
    #
    #         self.clustering_img = customtkinter.CTkImage(Image.open(f"{self.path_dict['temp']}/gene_cluster.png"),
    #                                                      size=(800, 320))
    #         self.clustering_img_l.configure(image=self.clustering_img)
    #
    # def run_chr_clustering(self):
    #     self.t4_progress_bar.set(0)
    #     i = 0
    #     for key, value in self.knowledge.items():
    #         self.t4_progress_bar.set((i + 1) / len(self.knowledge.items()))
    #         if self.file_name == value.chr_name.split("r")[1] and value.group == "cytobands":
    #             cgr = CGR(self.sequence[value.start:value.end], self.k_var.get())
    #             cgr_img = cgr.frequency_chaos_game_representation()
    #             cgr_img = normalize(cgr_img)
    #             self.cgr_images_clustering.append(cgr_img)
    #             self.band_info_clustering.append(value)
    #         i += 1


if __name__ == "__main__":
    app = App()
    app.mainloop()
