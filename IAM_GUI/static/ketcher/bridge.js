// Minimal bridge for postMessage round-trips during tests
(function () {
  window.addEventListener('message', (evt) => {
    try {
      if (!evt || !evt.data) return;
      const msg = evt.data;
      // Echo back recognizable responses for smoke tests
      if (msg.type === 'ping') {
        evt.source.postMessage({ type: 'pong' }, evt.origin);
      }
      // Add other lightweight handlers here if needed
    } catch (e) { /* no-op */ }
  });
})();
