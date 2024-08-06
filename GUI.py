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

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib import pyplot as plt

from sklearn import cluster
from sklearn.decomposition import PCA

from chaos_game_representation import CGR
from chromosomes_holder import ChromosomesHolder
from constants import DISTANCE_METRICS_LIST, SCIENTIFIC_NAMES, RESOLUTION_DICT
from distances.distance_metrics import get_dist

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
        self.seq = ""
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
        self.temp_output_path = "./.gui_temp_outputs"
        if not os.path.exists(self.temp_output_path):
            os.makedirs(self.temp_output_path)
        self.assets_path = "./assets"

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
                                                       fg_color="#707370", width=600, height=200)
        self.t1_display_frame.grid(row=1, column=1, padx=(5, 5), pady=(5, 5), sticky="nsew")

        # Designing the config frame (F1)
        for row in range(10):  # Increased row count to account for empty rows
            self.t1_config_frame.grid_rowconfigure(row, weight=1)

        # Creating frames for chromosome 1 and chromosome 2
        self.t1_chr_frame = customtkinter.CTkFrame(self.t1_config_frame, corner_radius=20)
        self.t1_chr_frame.grid(row=0, columnspan=2, padx=5, pady=5, sticky="ew")

        # Chromosomes Widget
        species_list = list(SCIENTIFIC_NAMES.keys())
        self.t1_ds = {'1': GUIDataStructure(), '2': GUIDataStructure()}
        self.t1_species_combobox = {}
        self.t1_chr_combobox = {}
        for i in range(2):
            label = customtkinter.CTkLabel(self.t1_chr_frame, text=f"Chromosome {i + 1}: ", font=self.header_font)
            label.grid(row=i * 3, column=0, sticky="w", padx=10, pady=(5, 0))
            specie_label = customtkinter.CTkLabel(self.t1_chr_frame, text="Species: ", font=self.header_font)
            specie_label.grid(row=(i * 3) + 1, column=0, sticky="w", padx=10)
            chr_label = customtkinter.CTkLabel(self.t1_chr_frame, text="Chromosome name: ", font=self.header_font)
            chr_label.grid(row=(i * 3) + 1, column=1, sticky="w", padx=10)
            self.t1_species_combobox[f"{i + 1}"] = customtkinter.CTkComboBox(self.t1_chr_frame, values=species_list,
                                                                             width=100,
                                                                             variable=self.t1_ds[str(i + 1)].specie,
                                                                             command=partial(self.specie_change_event,
                                                                                             f"{i + 1}"))

            self.t1_species_combobox[f"{i + 1}"].grid(row=(i * 3) + 2, column=0, sticky="w", padx=10, pady=(0, 10))
            self.t1_species_combobox[f"{i + 1}"].set("")

            self.t1_chr_combobox[f"{i + 1}"] = customtkinter.CTkComboBox(self.t1_chr_frame, values=[],
                                                                         variable=self.t1_ds[str(i + 1)].chromosome,
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

        # plot button
        self.t1_display_frame.grid_rowconfigure(0, weight=1)
        self.t1_display_frame.grid_columnconfigure(0, weight=1)
        t1_plot_button = customtkinter.CTkButton(self.t1_config_frame, text="Plot", width=120, command=self.t1_plot)
        t1_plot_button.grid(row=8, columnspan=2, pady=(10, 10))

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
                                          variable=self.t1_ds[str(i + 1)].annotation, state="disable")
            self.t1_parts_name_combobox[str(i + 1)].grid(row=(i * 2) + 1, column=0, padx=10, pady=_pad_size)

            seq_label_s = customtkinter.CTkLabel(self.t1_slider_frame, text='Start')
            seq_label_s.grid(row=(i * 2), column=1, padx=5, pady=_pad_size)
            seq_label_e = customtkinter.CTkLabel(self.t1_slider_frame, text='End')
            seq_label_e.grid(row=(i * 2) + 1, column=1, padx=5, pady=_pad_size)

            self.start_seq_scale[str(i + 1)] = ttk.Scale(self.t1_slider_frame, from_=0,
                                                         to=len(self.t1_ds[str(i + 1)].seq),
                                                         orient=customtkinter.HORIZONTAL, length=700,
                                                         variable=self.t1_ds[str(i + 1)].start_seq,
                                                         command=partial(self.sequence_value_change, str(i + 1)))
            self.start_seq_scale[str(i + 1)].grid(row=(i * 2), column=2, pady=_pad_size)

            self.t1_start_seq_entry[str(i + 1)] = customtkinter.CTkEntry(self.t1_slider_frame,
                                                                         textvariable=self.t1_ds[str(i + 1)].start_txt)
            self.t1_start_seq_entry[str(i + 1)].bind('<FocusOut>', partial(self.sequence_value_change, "3"))
            self.t1_start_seq_entry[str(i + 1)].bind('<Key-Return>', partial(self.sequence_value_change, "3"))
            self.t1_start_seq_entry[str(i + 1)].grid(row=(i * 2), column=3, padx=5, pady=_pad_size)
            seq_s_e_label = customtkinter.CTkLabel(self.t1_slider_frame, text='bp')
            seq_s_e_label.grid(row=(i * 2), column=4, pady=_pad_size)

            self.end_seq_scale[str(i + 1)] = ttk.Scale(self.t1_slider_frame, from_=0,
                                                       to=self.t1_ds[str(i + 1)].end_seq.get(),
                                                       orient=customtkinter.HORIZONTAL, length=700,
                                                       variable=self.t1_ds[str(i + 1)].end_seq,
                                                       command=partial(self.sequence_value_change, str(i + 1)))
            self.end_seq_scale[str(i + 1)].grid(row=(i * 2) + 1, column=2, pady=_pad_size)
            self.t1_end_seq_entry[str(i + 1)] = customtkinter.CTkEntry(self.t1_slider_frame,
                                                                       textvariable=self.t1_ds[str(i + 1)].end_txt)
            self.t1_end_seq_entry[str(i + 1)].bind('<FocusOut>', partial(self.sequence_value_change, "3"))
            self.t1_end_seq_entry[str(i + 1)].bind('<Key-Return>', partial(self.sequence_value_change, "3"))
            self.t1_end_seq_entry[str(i + 1)].grid(row=(i * 2) + 1, column=3, padx=5, pady=_pad_size)
            seq_s_e_label_d = customtkinter.CTkLabel(self.t1_slider_frame, text='bp')
            seq_s_e_label_d.grid(row=(i * 2) + 1, column=4, pady=_pad_size)

        '''
            Second Page (Consecutive Window Analyzer)
            Configuring 
        '''
        self.tabview.tab(tab_names[1]).grid_columnconfigure(0, weight=1)
        self.tabview.tab(tab_names[1]).grid_columnconfigure(1, weight=10)

        self.tabview.tab(tab_names[1]).grid_rowconfigure(0, weight=1)
        self.tabview.tab(tab_names[1]).grid_rowconfigure(1, weight=20)
        self.tabview.tab(tab_names[1]).grid_rowconfigure(2, weight=20)
        self.tabview.tab(tab_names[1]).grid_rowconfigure(3, weight=1)

        # Frames
        self.t2_config_frame = customtkinter.CTkFrame(self.tabview.tab(tab_names[1]), corner_radius=20,
                                                      fg_color="#333333")
        self.t2_config_frame.grid(row=0, column=0, rowspan=4, sticky="nsew")
        self.t2_plot_frame = customtkinter.CTkFrame(self.tabview.tab(tab_names[1]),
                                                    corner_radius=20, fg_color="white", width=1100, height=300)
        self.t2_plot_frame.grid(row=1, column=1, padx=(5, 5), pady=(5, 5), sticky="nsew")
        self.t2_display_frame = customtkinter.CTkFrame(self.tabview.tab(tab_names[1]), corner_radius=20,
                                                       fg_color="#707370", width=600, height=200)
        self.t2_display_frame.grid(row=2, column=1, padx=(5, 5), pady=(5, 5), sticky="nsew")

        self.t2_display_frame.grid_rowconfigure(0, weight=1)
        self.t2_display_frame.grid_columnconfigure(0, weight=1)

        self.t2_plot_frame.grid_rowconfigure(0, weight=1)
        self.t2_plot_frame.grid_columnconfigure(0, weight=1)

        # Designing the config frame (F2)
        for row in range(10):
            self.t2_config_frame.grid_rowconfigure(row, weight=1)

        # Creating frames for chromosome 1 and chromosome 2
        self.t2_chr_frame = customtkinter.CTkFrame(self.t2_config_frame, corner_radius=20)
        self.t2_chr_frame.grid(row=0, columnspan=2, padx=5, pady=5, sticky="ew")

        self.t2_ds = {'1': GUIDataStructure()}

        label = customtkinter.CTkLabel(self.t2_chr_frame, text=f"Chromosome: ", font=self.header_font)
        label.grid(row=0, column=0, sticky="w", padx=10, pady=(5, 0))
        specie_label = customtkinter.CTkLabel(self.t2_chr_frame, text="Species: ", font=self.header_font)
        specie_label.grid(row=1, column=0, sticky="w", padx=10)
        chr_label = customtkinter.CTkLabel(self.t2_chr_frame, text="Chromosome name: ", font=self.header_font)
        chr_label.grid(row=1, column=1, sticky="w", padx=10)
        self.t2_specie_combobox = customtkinter.CTkComboBox(self.t2_chr_frame, values=species_list, width=100,
                                                            variable=self.t2_ds["1"].specie,
                                                            command=partial(self.t2_specie_change_event))

        self.t2_specie_combobox.grid(row=2, column=0, sticky="w", padx=10, pady=(0, 10))
        self.t2_specie_combobox.set("")

        self.t2_chr_combobox = customtkinter.CTkComboBox(self.t2_chr_frame, values=[],
                                                         variable=self.t2_ds["1"].chromosome, width=100,
                                                         command=partial(self.t2_chromosome_change_event))
        self.t2_chr_combobox.grid(row=2, column=1, sticky="w", padx=10, pady=(0, 10))
        self.t2_chr_combobox.set("")

        # k_mer combo box
        k_mer_label = customtkinter.CTkLabel(self.t2_config_frame, text="k-mer: ", font=self.header_font)
        k_mer_label.grid(row=1, column=0, padx=10, pady=(10, 10))
        self.t2_k_mer_combobox = customtkinter.CTkComboBox(self.t2_config_frame, values=values_list, width=100,
                                                           state="normal", variable=self.k_var)
        self.t2_k_mer_combobox.grid(row=1, column=1, sticky="w", padx=10, pady=(10, 10))

        # Window size
        window_s_label = customtkinter.CTkLabel(self.t2_config_frame, text="Window Size:", font=self.header_font)
        window_s_label.grid(row=2, column=0, padx=10, pady=(20, 5))
        self.t2_window_s = tkinter.StringVar(value="")
        self.t2_window_entry = customtkinter.CTkEntry(self.t2_config_frame, textvariable=self.t2_window_s)
        self.t2_window_entry.configure(state="disable")
        self.t2_window_entry.grid(row=2, column=1, pady=(20, 5), padx=(0, 20), sticky="w")

        # Distance metrics
        t2_dist_metric_l = customtkinter.CTkLabel(self.t2_config_frame, text="Distance Metric: ", font=self.header_font)
        t2_dist_metric_l.grid(row=3, column=0, pady=(20, 5), padx=(5, 0))
        self.t2_dist_metric_combobox = customtkinter.CTkComboBox(self.t2_config_frame, values=DISTANCE_METRICS_LIST,
                                                                 width=120, variable=self.dist_metric)
        self.t2_dist_metric_combobox.grid(row=3, column=1, pady=(20, 5), sticky="w")
        self.t2_dist_metric_combobox.set("")

        t2_switch = customtkinter.CTkSwitch(self.t2_config_frame, text=f"Frequency CGR", variable=self.fcgr)
        t2_switch.grid(row=4, columnspan=2, pady=(20, 5))

        # run button
        t2_run_button = customtkinter.CTkButton(self.t2_config_frame, text="Run",
                                                command=partial(self.run_consecutive, None))
        t2_run_button.grid(row=7, columnspan=2)  # , sticky="ns")

        # cgr_img_base = customtkinter.CTkImage(Image.open(f"{self.assets_path}/background.png"), size=(800, 300))
        # self.t2_cgr_img_l = customtkinter.CTkLabel(self.t2_display_frame, image=cgr_img_base, text="")
        # self.t2_cgr_img_l.grid(row=0, padx=120, sticky='nsew')

        # placing the progress bars
        self.step_length = None
        self.cgr_distance_history = None
        self.progress_bar = customtkinter.CTkProgressBar(self.tabview.tab(tab_names[1]))
        self.progress_bar.set(0)
        self.progress_bar.grid(row=0, column=1, padx=(10, 10), pady=(10, 10), sticky="nsew")

        # placing the slider bar
        self.t2_changing_frame = customtkinter.CTkFrame(self.tabview.tab(tab_names[1]), corner_radius=20)
        self.t2_changing_frame.grid(row=3, column=1, sticky="nsew")

        # Designing the changing frame
        self.t2_changing_frame.grid_columnconfigure(0, weight=1)
        self.t2_changing_frame.grid_columnconfigure(1, weight=10)
        self.t2_changing_frame.grid_columnconfigure(2, weight=1)

        self.t2_pic_num = customtkinter.IntVar(value=0)
        self.t2_scale = customtkinter.CTkSlider(self.t2_changing_frame, from_=0,
                                                orientation=customtkinter.HORIZONTAL, variable=self.t2_pic_num,
                                                command=partial(self._change_images, self.t2_pic_num.get(), None))
        self.t2_scale.grid(row=0, column=1, pady=(10, 10), sticky="nsew")

        # previous-next button
        self.previous_im = customtkinter.CTkImage(light_image=Image.open(f"{self.assets_path}/back_arrow.png"),
                                                  size=(20, 20))
        self.next_im = customtkinter.CTkImage(light_image=Image.open(f"{self.assets_path}/next_arrow.png"),
                                              size=(20, 20))
        self.t2_previous_button = customtkinter.CTkButton(self.t2_changing_frame, image=self.previous_im, text="",
                                                          command=partial(self.move_previous, None), width=10)
        self.t2_previous_button.grid(row=0, column=0)

        self.t2_next_button = customtkinter.CTkButton(self.t2_changing_frame, image=self.next_im, text="",
                                                      command=partial(self.move_next, None), width=10)
        self.t2_next_button.grid(row=0, column=2)

    def sync_text_vars(self, sender, keep_annotation=False):
        self.t1_ds[sender].start_txt.set(f"{self.t1_ds[sender].start_seq.get()}")
        self.t1_ds[sender].end_txt.set(f"{self.t1_ds[sender].end_seq.get()}")
        if keep_annotation is False:
            self.t1_ds[sender].annotation.set("")

    def reverse_sync_text_vars(self, sender):
        self.t1_ds[sender].start_seq.set(int(self.t1_ds[sender].start_txt.get()))
        self.t1_ds[sender].end_seq.set(int(self.t1_ds[sender].end_txt.get()))
        self.t1_ds[sender].annotation.set("")

    # '''events'''
    def specie_change_event(self, sender, value):
        specie = self.t1_ds[sender].specie.get()
        self.t1_chr_combobox[sender].configure(values=ChromosomesHolder(specie).get_all_chromosomes_name())
        self.t1_ds[sender].invalidate_based_specie()

        # clear window_s
        self.window_s_toggle.set(0)
        self.window_size_toggle_event()

    def chromosome_change_event(self, sender, value):
        specie = self.t1_ds[sender].specie.get()
        chromosome = self.t1_ds[sender].chromosome.get()
        # set its sequence
        self.t1_ds[sender].seq = ChromosomesHolder(specie).get_chromosome_sequence(chromosome)
        # set the annotation combobox
        self.t1_parts_name_combobox[sender].configure(state="normal")
        self.t1_parts_name_combobox[sender].configure(values=ChromosomesHolder(specie).cytobands[chromosome])
        self.t1_ds[sender].invalidate_based_chromosome()
        # set start and end
        self.t1_ds[sender].end_seq.set(len(self.t1_ds[sender].seq))
        if len(self.t1_ds[sender].seq) > 0:
            self.start_seq_scale[sender].configure(state="normal")
            self.end_seq_scale[sender].configure(state="normal")
        else:
            self.start_seq_scale[sender].configure(state="disable")
            self.end_seq_scale[sender].configure(state="disable")
        self.start_seq_scale[sender].configure(to=len(self.t1_ds[sender].seq))
        self.end_seq_scale[sender].configure(to=len(self.t1_ds[sender].seq))

        self.sync_text_vars(sender)
        # clear window_s
        self.window_s_toggle.set(0)
        self.window_size_toggle_event()

    def annotation_change_event(self, sender, value):
        annotation = self.t1_ds[sender].annotation.get()
        chromosome_name = self.t1_ds[sender].chromosome.get()
        annotation_info = ChromosomesHolder(self.t1_ds[sender].specie.get()).cytobands[chromosome_name][annotation]
        self.t1_ds[sender].start_seq.set(annotation_info.start)
        self.t1_ds[sender].end_seq.set(annotation_info.end)
        self.sync_text_vars(sender, keep_annotation=True)

        # clear window_s
        self.window_s_toggle.set(0)
        self.window_size_toggle_event()

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
            for key, value in self.t1_ds.items():
                self.t1_ds[key].end_seq.set(self.t1_ds[key].start_seq.get() + int(self.window_s.get()))

        for key, value in self.t1_ds.items():
            self.sync_text_vars(key)

    def sequence_value_change(self, sender, value):
        if sender == "0":  # Window size changed
            for key, value in self.t1_ds.items():
                self.t1_ds[key].end_seq.set(self.t1_ds[key].start_seq.get() + int(self.window_s.get()))
        elif sender in ["1", "2"]:  # Scale changed
            if self.window_s_toggle.get() == 1:
                self.t1_ds[sender].end_seq.set(self.t1_ds[sender].start_seq.get() + int(self.window_s.get()))
        elif sender in ["3"]:
            for key, value in self.t1_ds.items():
                self.reverse_sync_text_vars(key)
            if self.window_s_toggle.get() == 1:
                for key, value in self.t1_ds.items():
                    self.t1_ds[key].end_seq.set(self.t1_ds[key].start_seq.get() + int(self.window_s.get()))

        for key, value in self.t1_ds.items():
            self.sync_text_vars(key)

    def t1_plot(self):
        fcgrs_dict = {}
        for key in self.t1_ds.keys():
            fcgrs_dict[key] = {}
            seq = self.t1_ds[key].seq[self.t1_ds[key].start_seq.get():self.t1_ds[key].end_seq.get()]
            if self.checkbox_RC[key].get():
                seq = ChromosomesHolder.get_reverse_complement(seq)
            if self.checkbox_Random[key].get():
                seq = list(seq)
                random.shuffle(seq)
                seq = ''.join(seq)
            cgr = CGR(seq, self.k_var.get())
            if self.fcgr.get() == 1:
                fcgrs_dict[key]["(f)cgr"] = cgr.get_fcgr()
            else:
                fcgrs_dict[key]["(f)cgr"] = cgr.get_cgr()

            fcgrs_dict[key]["chr_len"] = len(self.t1_ds[key].seq)
            fcgrs_dict[key]["b"] = self.t1_ds[key].start_seq.get()
            fcgrs_dict[key]["e"] = self.t1_ds[key].end_seq.get()

        diff = fcgrs_dict["2"]["(f)cgr"] - fcgrs_dict["1"]["(f)cgr"]
        fcgrs_dict["diff"] = diff
        distance_value = get_dist(fcgrs_dict["1"]["(f)cgr"], fcgrs_dict["2"]["(f)cgr"], dist_m=self.dist_metric.get())
        fcgrs_dict["distance"] = distance_value

        # Visualize the FCGRs
        display_frame_color = self.t1_display_frame.cget("fg_color")
        fig = self.plot_fcgrs(fcgrs_dict, colormap=True, background_color=display_frame_color)

        # Clear the previous figure from the display frame if any
        for widget in self.t1_display_frame.winfo_children():
            widget.destroy()

        # Create a canvas and add the figure to it
        canvas = FigureCanvasTkAgg(fig, master=self.t1_display_frame)
        canvas.draw()

        # Set the canvas size explicitly
        canvas_width = 600  # Example width, adjust as needed
        canvas_height = 200  # Example height, adjust as needed
        canvas.get_tk_widget().config(width=canvas_width, height=canvas_height)

        # Use grid to place the canvas
        canvas.get_tk_widget().grid(row=0, column=0, padx=10, pady=10, sticky='nsew')

    def plot_fcgrs(self, fcgrs, colormap=False, background_color=None):
        fig, (ax1, ax2, ax3) = plt.subplots(1, 3)
        extent = 0, 1, 0, 1

        if background_color is not None:
            fig.patch.set_facecolor(background_color)

        scale_1, scaling_1 = self.get_scaling(fcgrs["1"]["chr_len"])
        b1 = fcgrs["1"]["b"]
        e1 = fcgrs["1"]["e"]
        scale_2, scaling_2 = self.get_scaling(fcgrs["2"]["chr_len"])
        b2 = fcgrs["2"]["b"]
        e2 = fcgrs["2"]["e"]

        # plot the data on the subplots
        img1 = CGR.array2img(fcgrs["1"]["(f)cgr"], bits=8, resolution=RESOLUTION_DICT[self.k_var.get()])
        img1 = Image.fromarray(img1, 'L')
        ax1.imshow(img1, cmap='gray', extent=extent)
        ax1.tick_params(left=False, right=False, labelleft=False, labelbottom=False, bottom=False)
        ax1.set_title(f'{round(b1 / scale_1, 2)} - {round(e1 / scale_1, 2)} {scaling_1}')

        im2 = ax2.imshow(fcgrs['diff'], cmap='RdBu', norm=plt.Normalize(-100, 100), extent=extent)
        ax2.tick_params(left=False, right=False, labelleft=False, labelbottom=False, bottom=False)
        ax2.set_title(f'distance = {round(fcgrs["distance"], 4)}')

        img2 = CGR.array2img(fcgrs["2"]["(f)cgr"], bits=8, resolution=RESOLUTION_DICT[self.k_var.get()])
        img2 = Image.fromarray(img2, 'L')
        ax3.imshow(img2, cmap='gray', extent=extent)
        ax3.tick_params(left=False, right=False, labelleft=False, labelbottom=False, bottom=False)
        ax3.set_title(f'{round(b2 / scale_2, 2)} - {round(e2 / scale_2, 2)} {scaling_2}')

        if colormap:
            cbar_ax2 = fig.add_axes([0.36, 0.1, 0.3, 0.02])  # Adjust position as needed
            fig.colorbar(im2, cax=cbar_ax2, orientation='horizontal')

        return fig

    @staticmethod
    def get_scaling(chromosome_length):
        scale = 1_000_000
        while (chromosome_length / scale) < 2:
            scale //= 1000
        if scale == 1_000_000:
            scaling = "Mbp"
        elif scale == 1_000:
            scaling = "Kbp"
        else:
            scaling = "bp"
        return scale, scaling

    def t2_specie_change_event(self, value):
        specie = self.t2_ds["1"].specie.get()
        self.t2_chr_combobox.configure(values=ChromosomesHolder(specie).get_all_chromosomes_name())
        self.t2_ds["1"].invalidate_based_specie()

    def t2_chromosome_change_event(self, value):
        specie = self.t2_ds["1"].specie.get()
        chromosome = self.t2_ds["1"].chromosome.get()
        # set its sequence
        self.t2_ds["1"].seq = ChromosomesHolder(specie).get_chromosome_sequence(chromosome)
        # set end
        self.t2_ds["1"].end_seq.set(len(self.t2_ds["1"].seq))
        # enable window size
        self.t2_window_entry.configure(state="normal")

    def run_consecutive(self, event):
        global foo_thread
        foo_thread = threading.Thread(target=self.t2_run)
        foo_thread.daemon = True
        foo_thread.start()
        self.after(20, self.check_thread)

    def t2_run(self):
        self.cgr_distance_history = []
        self.step_length = np.floor(len(self.t2_ds["1"].seq) / int(self.t2_window_s.get()))
        self.progress_bar.set(0)
        self.t2_pic_num.set(0)
        for i in range(int(self.step_length) - 1):
            dictionary = {}
            self.progress_bar.set((i + 2) / int(self.step_length))
            b1 = i * int(self.t2_window_s.get())
            e1 = (i + 1) * int(self.t2_window_s.get())
            b2 = (i + 1) * int(self.t2_window_s.get())
            e2 = (i + 2) * int(self.t2_window_s.get())

            cgr1 = CGR(self.t2_ds["1"].seq[b1:e1], self.k_var.get())
            cgr2 = CGR(self.t2_ds["1"].seq[b2:e2], self.k_var.get())

            if self.fcgr.get() == 1:
                im1 = cgr1.get_fcgr()
                im2 = cgr2.get_fcgr()
            else:
                im1 = cgr1.get_cgr()
                im2 = cgr2.get_cgr()

            diff = im2 - im1

            dist = get_dist(im1, im2, dist_m=self.dist_metric.get())

            self.cgr_distance_history.append(dist)

            dictionary["1"] = {"(f)cgr": im1, "b": b1, "e": e1, "chr_len": len(self.t2_ds["1"].seq)}
            dictionary["2"] = {"(f)cgr": im2, "b": b2, "e": e2, "chr_len": len(self.t2_ds["1"].seq)}
            dictionary["diff"] = diff
            dictionary["distance"] = dist

            path = f"{self.temp_output_path}/consecutive/pickle"
            if not os.path.exists(path):
                os.makedirs(path)
            with open(f"{path}/{i}.pkl", 'wb') as f:
                pickle.dump(dictionary, f)

    def check_thread(self):
        if foo_thread.is_alive():
            self.after(20, self.check_thread)
        else:
            self.t2_scale.configure(to=int(len(self.cgr_distance_history) - 1))  # Update the scale range
            self._change_images(0, None)

    def _change_images(self, index, value):
        # plot distance results bar and first index is red
        self._t2_plot_chart(index)
        # Load and display the first image set in next plot
        self._t2_plot_fcgrs(index)

    def _t2_plot_chart(self, highlighted_index):
        fig, ax1 = plt.subplots(figsize=(50, 3))
        x = np.arange(len(self.cgr_distance_history))
        y = np.asarray(self.cgr_distance_history)
        mask1 = x == highlighted_index
        mask2 = x != highlighted_index
        # bar_width = 0.5
        ax1.bar(x[mask1], y[mask1], color='red')  # , width=bar_width)
        ax1.bar(x[mask2], y[mask2], color='blue')  # , width=bar_width)

        # Clear the previous figure from the display frame if any
        for widget in self.t2_display_frame.winfo_children():
            widget.destroy()

        # Create a canvas and add the figure to it
        canvas = FigureCanvasTkAgg(fig, master=self.t2_plot_frame)
        canvas.draw()

        # Set the canvas size explicitly
        canvas_width = 1100  # Example width, adjust as needed
        canvas_height = 300  # Example height, adjust as needed
        canvas.get_tk_widget().config(width=canvas_width, height=canvas_height)
        # Use grid to place the canvas
        canvas.get_tk_widget().grid(row=0, column=0, padx=(8, 8), pady=5, sticky='nsew')

    def _t2_plot_fcgrs(self, image_index, mode="consecutive"):
        with open(f"{self.temp_output_path}/{mode}/pickle/{image_index}.pkl", 'rb') as handle:
            dictionary = pickle.load(handle)

        display_frame_color = self.t2_display_frame.cget("fg_color")
        fig = self.plot_fcgrs(dictionary, colormap=False, background_color=display_frame_color)

        # Create a canvas and add the figure to it
        canvas = FigureCanvasTkAgg(fig, master=self.t2_display_frame)
        canvas.draw()

        # Set the canvas size explicitly
        canvas_width = 600  # Example width, adjust as needed
        canvas_height = 200  # Example height, adjust as needed
        canvas.get_tk_widget().config(width=canvas_width, height=canvas_height)
        # Use grid to place the canvas
        canvas.get_tk_widget().grid(row=0, column=0, padx=10, pady=10, sticky='nsew')

    def move_previous(self, value):
        if self.t2_pic_num.get() > 0:
            self.t2_pic_num.set(self.t2_pic_num.get() - 1)
            self._change_images(self.t2_pic_num.get(), None)
            # self._t2_plot_chart(self.t2_pic_num.get())
            # self._t2_plot_fcgrs(self.t2_pic_num.get())

    def move_next(self, value):
        if self.t2_pic_num.get() < len(self.cgr_distance_history) - 1:
            self.t2_pic_num.set(self.t2_pic_num.get() + 1)
            self._change_images(self.t2_pic_num.get(), None)
            # self._t2_plot_chart(self.t2_pic_num.get())
            # self._t2_plot_fcgrs(self.t2_pic_num.get())


if __name__ == "__main__":
    app = App()
    app.mainloop()
