// Offline cache. Bump VERSION when files change so phones pick up the update.
const VERSION = 'gyro-v6';
const FILES = ['./', 'index.html', 'physics.js', 'manifest.webmanifest',
  'icons/icon-180.png', 'icons/icon-192.png', 'icons/icon-512.png',
  'laser/', 'laser/index.html', 'laser/laser.js', 'laser/manifest.webmanifest',
  'laser/icons/icon-180.png', 'laser/icons/icon-192.png', 'laser/icons/icon-512.png',
  'qdots/', 'qdots/index.html', 'qdots/qd.js', 'qdots/manifest.webmanifest',
  'qdots/icons/icon-180.png', 'qdots/icons/icon-192.png', 'qdots/icons/icon-512.png'];
self.addEventListener('install', e => {
  e.waitUntil(caches.open(VERSION).then(c => c.addAll(FILES)).then(() => self.skipWaiting()));
});
self.addEventListener('activate', e => {
  e.waitUntil(caches.keys().then(ks => Promise.all(ks.filter(k => k !== VERSION).map(k => caches.delete(k))))
    .then(() => self.clients.claim()));
});
// network first, so updates arrive when online; cache when offline
self.addEventListener('fetch', e => {
  if (e.request.method !== 'GET') return;
  e.respondWith(fetch(e.request).then(r => {
    const copy = r.clone(); caches.open(VERSION).then(c => c.put(e.request, copy)); return r;
  }).catch(() => caches.match(e.request, { ignoreSearch: true })));
});
