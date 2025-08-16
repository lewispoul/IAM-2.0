// bridge.js - PostMessage bridge for Ketcher iframe communication

class KetcherBridge {
    constructor(iframeId) {
        this.iframe = document.getElementById(iframeId);
        this.ready = false;
        this.retryCount = 0;
        this.maxRetries = 20; // 10 seconds at 500ms intervals
        
        if (!this.iframe) {
            console.error('Ketcher iframe not found:', iframeId);
            return;
        }
        
        // Listen for Ketcher ready message
        window.addEventListener('message', (event) => {
            if (event.data && event.data.type === 'ketcher-ready') {
                this.ready = true;
                console.log('🧪 Ketcher bridge ready (via postMessage)');
                this.onReady();
            }
        });
        
        // Wait for iframe to load, then check for Ketcher
        this.iframe.addEventListener('load', () => {
            console.log('📄 Ketcher iframe loaded, checking for Ketcher...');
            this.checkKetcherReady();
        });
        
        // If iframe is already loaded, start checking immediately
        if (this.iframe.contentDocument && this.iframe.contentDocument.readyState === 'complete') {
            setTimeout(() => this.checkKetcherReady(), 100);
        }
    }
    
    checkKetcherReady() {
        try {
            const contentWindow = this.iframe.contentWindow;
            
            if (contentWindow && contentWindow.ketcher && typeof contentWindow.getMolfile === 'function') {
                this.ready = true;
                console.log('🧪 Ketcher bridge ready (direct access)');
                this.onReady();
                return;
            }
            
            // Retry if not ready yet
            if (this.retryCount < this.maxRetries) {
                this.retryCount++;
                setTimeout(() => this.checkKetcherReady(), 500);
            } else {
                console.warn('⚠️ Ketcher initialization timeout after', this.maxRetries * 0.5, 'seconds');
                this.onTimeout();
            }
            
        } catch (error) {
            // Cross-origin error - wait for postMessage
            if (this.retryCount < this.maxRetries) {
                this.retryCount++;
                setTimeout(() => this.checkKetcherReady(), 500);
            } else {
                console.warn('⚠️ Ketcher bridge timeout (cross-origin)');
                this.onTimeout();
            }
        }
    }
    
    onReady() {
        // Notify other parts of the application
        const event = new CustomEvent('ketcherReady', { detail: { bridge: this } });
        document.dispatchEvent(event);
    }
    
    onTimeout() {
        // Show timeout message
        const event = new CustomEvent('ketcherTimeout', { detail: { bridge: this } });
        document.dispatchEvent(event);
    }

    async getMolfile() {
        if (!this.ready) {
            throw new Error('Ketcher not ready. Please wait for initialization.');
        }
        
        try {
            // Direct access (same-origin)
            const contentWindow = this.iframe.contentWindow;
            if (contentWindow && contentWindow.getMolfile) {
                const molfile = contentWindow.getMolfile();
                return molfile || '';
            }
            
            // Fallback: try accessing ketcher directly
            if (contentWindow && contentWindow.ketcher && contentWindow.ketcher.getMolfile) {
                const molfile = await contentWindow.ketcher.getMolfile();
                return molfile || '';
            }
            
            throw new Error('getMolfile function not available');
            
        } catch (error) {
            console.error('Failed to get molfile:', error);
            throw new Error(`Could not get molfile: ${error.message}`);
        }
    }

    async setMolfile(molfile) {
        if (!this.ready) {
            throw new Error('Ketcher not ready. Please wait for initialization.');
        }
        
        try {
            const contentWindow = this.iframe.contentWindow;
            if (contentWindow && contentWindow.setMolfile) {
                return contentWindow.setMolfile(molfile);
            }
            
            if (contentWindow && contentWindow.ketcher && contentWindow.ketcher.setMolfile) {
                return await contentWindow.ketcher.setMolfile(molfile);
            }
            
            throw new Error('setMolfile function not available');
            
        } catch (error) {
            console.error('Failed to set molfile:', error);
            throw new Error(`Could not set molfile: ${error.message}`);
        }
    }

    async getSmiles() {
        if (!this.ready) {
            throw new Error('Ketcher not ready. Please wait for initialization.');
        }
        
        try {
            const contentWindow = this.iframe.contentWindow;
            if (contentWindow && contentWindow.getSmiles) {
                const smiles = contentWindow.getSmiles();
                return smiles || '';
            }
            
            if (contentWindow && contentWindow.ketcher && contentWindow.ketcher.getSmiles) {
                const smiles = await contentWindow.ketcher.getSmiles();
                return smiles || '';
            }
            
            throw new Error('getSmiles function not available');
            
        } catch (error) {
            console.error('Failed to get SMILES:', error);
            throw new Error(`Could not get SMILES: ${error.message}`);
        }
    }

    async setSmiles(smiles) {
        if (!this.ready) {
            throw new Error('Ketcher not ready. Please wait for initialization.');
        }
        
        try {
            const contentWindow = this.iframe.contentWindow;
            if (contentWindow && contentWindow.setSmiles) {
                return contentWindow.setSmiles(smiles);
            }
            
            if (contentWindow && contentWindow.ketcher && contentWindow.ketcher.setSmiles) {
                return await contentWindow.ketcher.setSmiles(smiles);
            }
            
            throw new Error('setSmiles function not available');
            
        } catch (error) {
            console.error('Failed to set SMILES:', error);
            throw new Error(`Could not set SMILES: ${error.message}`);
        }
    }

    async clear() {
        if (!this.ready) {
            throw new Error('Ketcher not ready. Please wait for initialization.');
        }
        
        try {
            const contentWindow = this.iframe.contentWindow;
            if (contentWindow && contentWindow.clear) {
                return contentWindow.clear();
            }
            
            if (contentWindow && contentWindow.ketcher) {
                return await contentWindow.ketcher.setMolfile('');
            }
            
            throw new Error('clear function not available');
            
        } catch (error) {
            console.error('Failed to clear:', error);
            throw new Error(`Could not clear editor: ${error.message}`);
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
        console.log('🔧 Ketcher bridge initialized');
    } else {
        console.error('❌ Ketcher iframe not found');
    }
});

// Export for debugging
window.KetcherBridge = KetcherBridge;
