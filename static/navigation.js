(() => {
  const header = document.querySelector('.app-header');
  const toggle = header.querySelector('.nav-toggle');
  const navigation = header.querySelector('.header-nav');
  const mobile = window.matchMedia('(max-width: 760px)');

  const setOpen = (open) => {
    header.classList.toggle('nav-open', open);
    toggle.setAttribute('aria-expanded', String(open));
    toggle.setAttribute('aria-label', open ? 'Close navigation' : 'Open navigation');
  };

  toggle.addEventListener('click', () => {
    setOpen(toggle.getAttribute('aria-expanded') !== 'true');
  });

  navigation.addEventListener('click', (event) => {
    if (mobile.matches && event.target.closest('a')) {
      setOpen(false);
      toggle.focus();
    }
  });

  document.addEventListener('click', (event) => {
    if (!header.contains(event.target)) setOpen(false);
  });

  document.addEventListener('keydown', (event) => {
    if (event.key === 'Escape' && toggle.getAttribute('aria-expanded') === 'true') {
      setOpen(false);
      toggle.focus();
    }
  });

  mobile.addEventListener('change', () => {
    const focused = document.activeElement;
    setOpen(false);
    if (mobile.matches && navigation.contains(focused)) toggle.focus();
    if (!mobile.matches && focused === toggle) navigation.querySelector('a').focus();
  });

  header.classList.add('nav-ready');
})();
