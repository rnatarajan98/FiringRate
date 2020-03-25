from dataclasses import dataclass
import matlab as matlab

@dataclass
class cell:
    name: str = None
    Spont_mean: float = None
    Spont_SEM: float = None
    rf_sta_lat: float = None
    rf_xdeg: float = None
    rf_ydeg: float = None
    rf_center: float = None
    rf_sigma: float = None
    sf_tuning: float = None
    sf_preferred: float = None
    sf_cutoff: float = None
    sf_latency: float = None
    sf_spikes: float = None
    tf_tuning: float = None
    tf_spikes: float = None
    c_tuning: float = None
    c_gain: float = None
    c_c50: float = None
    c_fit: float = None
    c_spikes: float = None
    d_tuning: float = None
    d_preferred: float = None
    d_DSI: float = None
    d_OSI: float = None
    d_spikes: float = None
    f_transient_infex: float = None
    f_spikes: float = None
    sf_curvefit: float = None
    c_curvefit: float = None
    
    
    