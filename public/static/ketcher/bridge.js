// bridge.js - PostMessage bridge for Ketcher iframe communication

class KetcherBridge {
    constructor(iframeId) {
        this.iframe = document.getElementById(iframeId);
        this.ready = false;
        
        // Wait for iframe to load
        this.iframe.addEventListener('load', () => {
            this.ready = true;
            console.log('🧪 Ketcher bridge ready');
        });
    }

    async getMolfile() {
        if (!this.ready || !this.iframe.contentWindow) {
            throw new Error('Ketcher not ready');
        }
        
        try {
            // Direct access since we're same-origin
            const molfile = this.iframe.contentWindow.getMolfile();
            return molfile || '';
        } catch (error) {
            console.error('Failed to get molfile:', error);
            throw error;
        }
    }

    async setMolfile(molfile) {
        if (!this.ready || !this.iframe.contentWindow) {
            throw new Error('Ketcher not ready');
        }
        
        try {
            return this.iframe.contentWindow.setMolfile(molfile);
        } catch (error) {
            console.error('Failed to set molfile:', error);
            throw error;
        }
    }

    async getSmiles() {
        if (!this.ready || !this.iframe.contentWindow) {
            throw new Error('Ketcher not ready');
        }
        
        try {
            const smiles = this.iframe.contentWindow.getSmiles();
            return smiles || '';
        } catch (error) {
            console.error('Failed to get SMILES:', error);
            throw error;
        }
    }

    async setSmiles(smiles) {
        if (!this.ready || !this.iframe.contentWindow) {
            throw new Error('Ketcher not ready');
        }
        
        try {
            return this.iframe.contentWindow.setSmiles(smiles);
        } catch (error) {
            console.error('Failed to set SMILES:', error);
            throw error;
        }
    }

    async clear() {
        if (!this.ready || !this.iframe.contentWindow) {
            throw new Error('Ketcher not ready');
        }
        
        try {
            return this.iframe.contentWindow.clear();
        } catch (error) {
            console.error('Failed to clear:', error);
            throw error;
        }
    }
}

// Global bridge instance
window.ketcherBridge = null;

// Initialize bridge when DOM is ready
document.addEventListener('DOMContentLoaded', () => {
    const iframe = document.getElementById('ketcherFrame');
    if (iframe) {
        window.ketcherBridge = new KetcherBridge('ketcherFrame');
    }
});
