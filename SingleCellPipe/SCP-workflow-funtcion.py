def RunSCVELO(adata=None, h5ad=None, group_by=None, neighbors_reduction=None, compare_reduction=None, dirpath="",fileprefix="",dpi=300,
              min_shared_counts=30, n_pcs=30, n_neighbors=30, approx=True, stream_smooth=0.3, stream_density=1.2,
              arrow_density=0.05, arrow_length=15, arrow_size=15, paga_threshold=0.15,
              calculate_velocity_genes=False, velocity_genes_min_corr=0.3, velocity_ngenes=100,
              s_genes=None, g2m_genes=None,
              recover_dynamics=False, n_jobs=12, velocity_with_noise=False,
              calculate_dynamical_genes=False, dynamical_ngenes=100,
              diff_kinetics=False):
  # import matplotlib
  # matplotlib.use('agg')
  import matplotlib.pyplot as plt
  import random
  random.seed(11)
  import scvelo as scv
  import pandas as pd
  import os
  prevdir = os.getcwd()
  os.chdir(os.path.expanduser(dirpath))

  try:
    if adata is None and h5ad is None:
      print("adata or h5ad must be provided.")
      exit()
    if group_by is None or neighbors_reduction is None or compare_reduction is None:
      print("group_by, neighbors_reduction and compare_reduction must be all provided.")
      exit()
    if adata is None:
      adata = scv.read(h5ad)
    del adata.uns
  
    neighbors_reduction = "X_" + neighbors_reduction
    compare_reduction = "X_" + compare_reduction
    adata.obs[group_by] = adata.obs[group_by].astype(dtype="category")
  
    scv.pp.filter_and_normalize(adata, min_shared_counts=min_shared_counts)
    scv.pp.moments(adata, n_pcs=n_pcs, use_rep=neighbors_reduction,
                   n_neighbors=n_neighbors)
  
    scv.tl.velocity(adata, vkey="stochastic")
    scv.tl.velocity_graph(adata, vkey="stochastic",
                          n_neighbors=n_neighbors, approx=approx)
    scv.pl.velocity_embedding_stream(adata, title="stochastic", basis=compare_reduction, vkey=[
                                     "stochastic"], color=group_by, smooth=stream_smooth, density=stream_density,save=False, show=True)
    plt.savefig('.'.join(filter(None,[fileprefix,"stochastic_stream.png"])),dpi=dpi)
    scv.pl.velocity_embedding(adata, title="stochastic", basis=compare_reduction, vkey=[
                              "stochastic"], color=group_by, arrow_length=arrow_length, arrow_size=arrow_size, density=arrow_density,save=False, show=True)
    plt.savefig('.'.join(filter(None,[fileprefix,"stochastic_arrow.png"])),dpi=dpi)
    
    scv.tl.velocity_confidence(adata, vkey="stochastic")
    scv.tl.velocity_pseudotime(adata, vkey="stochastic")
    scv.pl.scatter(adata, basis=compare_reduction, color=('stochastic_length', 'stochastic_confidence'),
                   cmap='coolwarm', perc=[5, 95], save=False, show=True)
    plt.savefig('.'.join(filter(None,[fileprefix,"stochastic_length_confidence.png"])),dpi=dpi)
    scv.pl.scatter(adata, basis=compare_reduction, color='stochastic_pseudotime',
                   cmap='gnuplot', save=False, show=True)
    plt.savefig('.'.join(filter(None,[fileprefix,"stochastic_pseudotime.png"])),dpi=dpi)
  
    adata.uns['neighbors']['distances'] = adata.obsp['distances']
    adata.uns['neighbors']['connectivities'] = adata.obsp['connectivities']
    scv.tl.paga(adata, groups=group_by, vkey="stochastic")
    scv.pl.paga(adata, basis=compare_reduction[2:], threshold=paga_threshold, size=50, alpha=0.02,
                min_edge_width=2, node_size_scale=1.5, save=False, show=True)
    plt.savefig('.'.join(filter(None,[fileprefix,"stochastic_paga.png"])),dpi=dpi)
  
    if calculate_velocity_genes is True:
      scv.tl.rank_velocity_genes(adata, vkey="stochastic", groupby=group_by,
                                 min_corr=velocity_genes_min_corr, n_genes=velocity_ngenes)
      df = scv.DataFrame(adata.uns['rank_velocity_genes']['names'])
      for cluster in df.columns:
        #df[0:1].values.ravel()[:12] ### by row
        scv.pl.scatter(adata, color=group_by, basis=df[cluster].values[:6], size=20, linewidth=2, alpha=1, ylabel="cluster: "+cluster+"\nunspliced",
                       add_linfit=True, add_rug=True, add_outline=True, ncols=3, frameon=True, save=False,show=False)
        plt.savefig('.'.join(filter(None,[fileprefix,cluster,"stochastic_genes1.png"])),dpi=dpi)
        scv.pl.velocity(adata, color=group_by, var_names=df[cluster].values[:6], size=10, linewidth=2, alpha=1, ylabel="cluster: "+cluster+"\nunspliced",
                        add_outline=True, basis=compare_reduction, color_map=["Spectral", "YlOrRd"], ncols=2,save=False,show=False)
        plt.savefig('.'.join(filter(None,[fileprefix,cluster,"stochastic_genes2.png"])),dpi=dpi)
        
    if s_genes is not None and g2m_genes is not None:
      scv.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)
      scv.pl.scatter(adata, basis=compare_reduction, color=('S_score', 'G2M_score'), save=False, show=True)
      plt.savefig('.'.join(filter(None,[fileprefix,"stochastic_cellcycle.png"])),dpi=dpi)
  
    if recover_dynamics is True or diff_kinetics is True or velocity_with_noise is True:
      adata2 = adata[:, adata.var['stochastic_genes']].copy()
      Ms = adata2.layers["Ms"]
      Mu = adata2.layers["Mu"]
      spliced = adata2.layers["spliced"]
      unspliced = adata2.layers["unspliced"]
      stochastic = adata2.layers["stochastic"]
      variance_stochastic = adata2.layers["variance_stochastic"]
      adata2.layers.clear()
      adata2.layers["Ms"] = Ms
      adata2.layers["Mu"] = Mu
      connectivities = adata2.obsp["connectivities"]
      distances = adata2.obsp["distances"]
      adata2.obsp.clear()
      adata2.obsp["connectivities"] = connectivities
  
      scv.tl.recover_dynamics(
          adata2, var_names='stochastic_genes', use_raw=False, n_jobs=n_jobs)
      adata2.obsp["distances"] = distances
      adata2.layers["spliced"] = spliced
      adata2.layers["unspliced"] = unspliced
      adata2.layers["stochastic"] = stochastic
      adata2.layers["variance_stochastic"] = variance_stochastic
      scv.tl.velocity(adata2, mode="dynamical", vkey="dynamical")
      scv.tl.velocity_graph(adata2, vkey="dynamical",
                            n_neighbors=n_neighbors, approx=approx)
      scv.pl.velocity_embedding_stream(adata2, title="dynamical", basis=compare_reduction, vkey=[
                                       "dynamical"], color=group_by, smooth=stream_smooth, density=stream_density, save=False, show=True)
      plt.savefig('.'.join(filter(None,[fileprefix,"dynamical_stream.png"])),dpi=dpi)
      scv.pl.velocity_embedding(adata2, title="dynamical", basis=compare_reduction, vkey=[
                                "dynamical"], color=group_by, arrow_length=arrow_length, arrow_size=arrow_size, density=arrow_density,save=False, show=True)
      plt.savefig('.'.join(filter(None,[fileprefix,"dynamical_arrow.png"])),dpi=dpi)
  
      scv.tl.velocity_confidence(adata2, vkey="dynamical")
      scv.tl.velocity_pseudotime(adata2, vkey="dynamical")
      scv.pl.scatter(adata2, basis=compare_reduction, color=('dynamical_length', 'dynamical_confidence'),
                     cmap='coolwarm', perc=[5, 95],save=False, show=True)
      plt.savefig('.'.join(filter(None,[fileprefix,"dynamical_length_confidence.png"])),dpi=dpi)
      scv.pl.scatter(adata2, basis=compare_reduction, color='dynamical_pseudotime',
                     cmap='gnuplot', save=False, show=True)
      plt.savefig('.'.join(filter(None,[fileprefix,"dynamical_pseudotime.png"])),dpi=dpi)
  
      scv.tl.latent_time(adata2, vkey="dynamical")
      scv.pl.scatter(adata2, basis=compare_reduction, color='latent_time',
                     color_map='gnuplot', save=False, show=True)
      plt.savefig('.'.join(filter(None,[fileprefix,"dynamical_latent_time.png"])),dpi=dpi)
  
      if calculate_dynamical_genes is True:
        scv.tl.rank_dynamical_genes(
            adata2, groupby=group_by, n_genes=dynamical_ngenes)
        df = scv.DataFrame(adata2.uns['rank_dynamical_genes']['names'])
        for cluster in df.columns:
          #df[0:1].values.ravel()[:12] ### by row
          scv.pl.scatter(adata, color=group_by, basis=df[cluster].values[:6], size=20, linewidth=2, alpha=1, ylabel="cluster: "+cluster+"\nunspliced",
                         add_linfit=True, add_rug=True, add_outline=True, ncols=3, frameon=True,save=False,show=False)
          plt.savefig('.'.join(filter(None,[fileprefix,cluster,"dynamical_genes1.png"])),dpi=dpi)
          scv.pl.velocity(adata, color=group_by, var_names=df[cluster].values[:6], size=10, linewidth=2, alpha=1, ylabel="cluster: "+cluster+"\nunspliced",
                          add_outline=True, basis=compare_reduction, color_map=["Spectral", "YlOrRd"], ncols=2, save=False,show=False)
          plt.savefig('.'.join(filter(None,[fileprefix,cluster,"dynamical_genes2.png"])),dpi=dpi)
          
      if diff_kinetics is True:
        top_genes = adata2.var['fit_likelihood'].sort_values(
            ascending=False).index[:100]
        scv.tl.differential_kinetic_test(
            adata2, var_names=top_genes, groupby=group_by)
        scv.tl.velocity(adata2, mode="dynamical",
                        vkey="dynamical_kinetics", diff_kinetics=True)
        scv.tl.velocity_graph(adata2, vkey="dynamical_kinetics",
                              n_neighbors=n_neighbors, approx=approx)
        scv.pl.velocity_embedding_stream(adata2, title="dynamical_kinetics", basis=compare_reduction, vkey=[
                                         "dynamical_kinetics"], color=group_by, smooth=stream_smooth, density=stream_density, save=False, show=True)
        plt.savefig('.'.join(filter(None,[fileprefix,"dynamical_kinetics_stream.png"])),dpi=dpi)
        scv.pl.velocity_embedding(adata2, title="dynamical_kinetics", basis=compare_reduction, vkey=[
                                  "dynamical_kinetics"], color=group_by, arrow_length=arrow_length, arrow_size=arrow_size, density=arrow_density, save=False, show=True)
        plt.savefig('.'.join(filter(None,[fileprefix,"dynamical_kinetics_arrow.png"])),dpi=dpi)
        scv.tl.velocity(adata2, mode="stochastic",
                        vkey="stochastic_kinetics", diff_kinetics=True)
        scv.tl.velocity_graph(adata2, vkey="stochastic_kinetics",
                              n_neighbors=n_neighbors, approx=approx)
        scv.pl.velocity_embedding_stream(adata2, title="stochastic_kinetics", basis=compare_reduction, vkey=[
                                         "stochastic_kinetics"], color=group_by, smooth=stream_smooth, density=stream_density, save=False, show=True)
        plt.savefig('.'.join(filter(None,[fileprefix,"stochastic_kinetics_stream.png"])),dpi=dpi)
        scv.pl.velocity_embedding(adata2, title="stochastic_kinetics", basis=compare_reduction, vkey=[
                                  "stochastic_kinetics"], color=group_by, arrow_length=arrow_length, arrow_size=arrow_size, density=arrow_density, save=False, show=True)
        plt.savefig('.'.join(filter(None,[fileprefix,"stochastic_kinetics_arrow.png"])),dpi=dpi)
        
      if velocity_with_noise is True:
        import numpy as np
        top_genes = adata2.var['fit_likelihood'].sort_values(
            ascending=False).index[:3]
        adata2.layers['dynamical_with_noise'] = adata2.layers['dynamical'] + \
            np.random.normal(
            adata2.layers['dynamical'], scale=adata2.layers['Ms'].std(0))
        scv.tl.velocity_graph(adata2, gene_subset=top_genes,
                              vkey='dynamical_with_noise')
        scv.tl.velocity_embedding(
            adata2, basis=compare_reduction[2:], vkey='dynamical_with_noise', autoscale=False)
        scv.pl.velocity_embedding_stream(adata2, title="dynamical_with_noise", basis=compare_reduction, vkey=[
                                         "dynamical_with_noise"], color=group_by, smooth=stream_smooth, density=stream_density,save=False, show=True)
        plt.savefig('.'.join(filter(None,[fileprefix,"dynamical_with_noise_stream.png"])),dpi=dpi)
        scv.pl.velocity_embedding(adata2, title="dynamical_with_noise", basis=compare_reduction, vkey=[
                                  "dynamical_with_noise"], color=group_by, arrow_length=arrow_length, arrow_size=arrow_size, density=arrow_density,save=False, show=True)
        plt.savefig('.'.join(filter(None,[fileprefix,"dynamical_with_noise_arrow.png"])),dpi=dpi)
        
      import os
      os.makedirs("data", exist_ok=True)
      adata2.write("data/"+prefix+".dynamical.h5ad", compression='gzip')
  finally:
        os.chdir(prevdir)
  
  try:
    adata.__dict__['_raw'].__dict__['_var'] = adata.__dict__[
        '_raw'].__dict__['_var'].rename(columns={'_index': 'features'})
  except:
    pass

  return adata


def RunPAGA(adata=None, h5ad=None, group_by=None, neighbors_reduction=None, compare_reduction=None, dirpath="",fileprefix="",dpi=300,
            n_pcs=30, threshold=0.15, use_rna_velocity=False):
  # import matplotlib
  # matplotlib.use('agg')
  import matplotlib.pyplot as plt
  import random
  random.seed(11)
  import scanpy as sc
  import os
  prevdir = os.getcwd()
  os.chdir(os.path.expanduser(dirpath))
  
  try:
    if adata is None and h5ad is None:
      print("adata or h5ad must be provided.")
      exit()
    if group_by is None or neighbors_reduction is None or compare_reduction is None:
      print("group_by, neighbors_reduction and compare_reduction must be all provided.")
      exit()
    if adata is None:
      adata = sc.read(h5ad)
    del adata.uns
  
    neighbors_reduction = "X_" + neighbors_reduction
    compare_reduction = "X_" + compare_reduction
    adata.obs[group_by] = adata.obs[group_by].astype(dtype="category")
  
    sc.pp.neighbors(adata, n_pcs=n_pcs, use_rep=neighbors_reduction)
    sc.tl.paga(adata, groups=group_by, use_rna_velocity=use_rna_velocity)
    if use_rna_velocity is True:
      sc.pl.paga(adata, threshold=threshold, arrowsize=10, transitions="transitions_confidence",
                 dashed_edges="connectivities", edge_width_scale=0.5, save=False, plot=False)
      sc.pl.paga_compare(adata, basis=compare_reduction, threshold=threshold, transitions="transitions_confidence",
                         dashed_edges="connectivities", edge_width_scale=0.5, save=False, show=True)
    else:
      sc.pl.paga(adata, threshold=threshold,
                 edge_width_scale=0.5, save=False, plot=False)
      sc.pl.paga_compare(adata, basis=compare_reduction, threshold=threshold,
                         edge_width_scale=0.5, save=False, show=True)
    plt.savefig('.'.join(filter(None,[fileprefix,"paga_compare.png"])),dpi=dpi)
  
    # sc.tl.leiden(adata, resolution = 1.0)
    # sc.tl.paga(adata, groups = "leiden")
    # sc.pl.paga(adata)
    # sc.tl.umap(adata, init_pos = "paga")
    # sc.pl.paga_compare(adata,basis = "umap",threshold = 0.15,edge_width_scale = 0.5,save="_tmp.png")
  finally:
        os.chdir(prevdir)
  
  try:
    adata.__dict__['_raw'].__dict__['_var'] = adata.__dict__[
        '_raw'].__dict__['_var'].rename(columns={'_index': 'features'})
  except:
    pass
  
  return adata
