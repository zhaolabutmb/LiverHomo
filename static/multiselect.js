document.addEventListener('DOMContentLoaded', () => {
  document.querySelectorAll('select[multiple]').forEach((select) => {
    const selectedValues = () => new Set(Array.from(select.selectedOptions, (option) => option.value));
    let previousValues = selectedValues();
    const rememberSelection = () => {
      previousValues = selectedValues();
    };

    select.addEventListener('focus', rememberSelection);
    select.addEventListener('change', () => {
      const allOption = Array.from(select.options).find((option) => option.value === '__all__');
      if (allOption?.selected && select.selectedOptions.length > 1) {
        if (previousValues.has(allOption.value)) {
          allOption.selected = false;
        } else {
          Array.from(select.options).forEach((option) => {
            option.selected = option === allOption;
          });
        }
      }
      rememberSelection();
    }, { capture: true });

    const toggleMouseOption = (event) => {
      const option = event.target instanceof HTMLOptionElement ? event.target : null;
      if (!option || option.matches(':disabled') || event.button !== 0) return;

      event.preventDefault();
      select.focus();

      const allOption = Array.from(select.options).find((item) => item.value === '__all__');
      if (option === allOption) {
        const selectAll = !option.selected;
        Array.from(select.options).forEach((item) => {
          item.selected = false;
        });
        option.selected = selectAll;
      } else {
        option.selected = !option.selected;
        if (allOption) allOption.selected = false;
      }

      select.dispatchEvent(new Event('change', { bubbles: true }));
    };

    if ('PointerEvent' in window) {
      select.addEventListener('pointerdown', (event) => {
        rememberSelection();
        if (event.pointerType === 'mouse') toggleMouseOption(event);
      });
    } else {
      // Older touch browsers must retain native selection rather than receive mouse emulation.
      select.addEventListener('mousedown', (event) => {
        rememberSelection();
        if (window.matchMedia('(hover: hover) and (pointer: fine)').matches) toggleMouseOption(event);
      });
    }
  });
});
