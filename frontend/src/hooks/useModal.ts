/**
 * useModal Hook
 * Reusable hook for modal dialogs with form integration
 * Eliminates code duplication across components with create/edit modals
 */

import { useState, useCallback } from 'react'
import type { FormInstance } from 'antd'
import { Form, message } from 'antd'

/**
 * Configuration options for useModal hook
 */
export interface UseModalOptions<T> {
  /**
   * Callback when form is submitted
   * Should handle create or update based on whether editing item exists
   */
  onSubmit?: (values: T, editing?: T | null) => Promise<void>

  /**
   * Callback when modal is opened
   */
  onOpen?: (item?: T) => void

  /**
   * Callback when modal is closed
   */
  onClose?: () => void

  /**
   * Callback when form submission succeeds
   */
  onSuccess?: () => void

  /**
   * Callback when form submission fails
   */
  onError?: (error: Error) => void

  /**
   * Whether to show success message after submission
   * @default true
   */
  showSuccessMessage?: boolean

  /**
   * Whether to show error message on submission failure
   * @default true
   */
  showErrorMessage?: boolean

  /**
   * Custom success messages
   */
  successMessages?: {
    create?: string
    update?: string
  }

  /**
   * Custom error messages
   */
  errorMessages?: {
    create?: string
    update?: string
  }

  /**
   * Whether to reset form on close
   * @default true
   */
  resetOnClose?: boolean
}

/**
 * Return type for useModal hook
 */
export interface UseModalReturn<T> {
  // State
  visible: boolean
  editing: T | null
  submitting: boolean
  form: FormInstance<T>

  // Actions
  open: (item?: T) => void
  close: () => void
  submit: (values: T) => Promise<void>
  handleOk: () => void
  handleCancel: () => void
}

/**
 * Hook for modal dialogs with form integration
 * Handles modal visibility, form state, and submission
 *
 * @example
 * const modal = useModal<ProjectCreate>({
 *   onSubmit: async (values, editing) => {
 *     if (editing) {
 *       await projectService.update(editing.id, values)
 *     } else {
 *       await projectService.create(values)
 *     }
 *   },
 *   onSuccess: () => refresh()
 * })
 *
 * // In component
 * <Button onClick={() => modal.open()}>Create</Button>
 * <Button onClick={() => modal.open(project)}>Edit</Button>
 *
 * <Modal
 *   open={modal.visible}
 *   onOk={modal.handleOk}
 *   onCancel={modal.handleCancel}
 *   confirmLoading={modal.submitting}
 * >
 *   <Form form={modal.form}>...</Form>
 * </Modal>
 */
export function useModal<T = any>(options: UseModalOptions<T> = {}): UseModalReturn<T> {
  const {
    onSubmit,
    onOpen,
    onClose,
    onSuccess,
    onError,
    showSuccessMessage = true,
    showErrorMessage = true,
    successMessages = {},
    errorMessages = {},
    resetOnClose = true,
  } = options

  // State
  const [visible, setVisible] = useState(false)
  const [editing, setEditing] = useState<T | null>(null)
  const [submitting, setSubmitting] = useState(false)
  const [form] = Form.useForm<T>()

  /**
   * Open modal (for create or edit)
   */
  const open = useCallback(
    (item?: T) => {
      if (item) {
        // Edit mode
        setEditing(item)
        form.setFieldsValue(item)
      } else {
        // Create mode
        setEditing(null)
        form.resetFields()
      }

      setVisible(true)
      onOpen?.(item)
    },
    [form, onOpen],
  )

  /**
   * Close modal
   */
  const close = useCallback(() => {
    setVisible(false)
    setEditing(null)

    if (resetOnClose) {
      form.resetFields()
    }

    onClose?.()
  }, [form, resetOnClose, onClose])

  /**
   * Submit form
   */
  const submit = useCallback(
    async (values: T) => {
      if (!onSubmit) {
        console.warn('onSubmit callback not provided')
        close()
        return
      }

      setSubmitting(true)

      try {
        await onSubmit(values, editing)

        // Show success message
        if (showSuccessMessage) {
          const defaultMessage = editing ? 'Updated successfully' : 'Created successfully'
          const customMessage = editing ? successMessages.update : successMessages.create

          message.success(customMessage || defaultMessage)
        }

        onSuccess?.()
        close()
      } catch (error: any) {
        // Show error message
        if (showErrorMessage) {
          const defaultMessage = editing ? `Failed to update: ${error.message}` : `Failed to create: ${error.message}`

          const customMessage = editing ? errorMessages.update : errorMessages.create

          message.error(customMessage || defaultMessage)
        }

        onError?.(error)
      } finally {
        setSubmitting(false)
      }
    },
    [
      onSubmit,
      editing,
      close,
      showSuccessMessage,
      showErrorMessage,
      successMessages,
      errorMessages,
      onSuccess,
      onError,
    ],
  )

  /**
   * Handle modal OK button click
   * Validates and submits form
   */
  const handleOk = useCallback(() => {
    form
      .validateFields()
      .then((values) => {
        submit(values)
      })
      .catch((info) => {
        console.log('Validation failed:', info)
      })
  }, [form, submit])

  /**
   * Handle modal cancel button click
   */
  const handleCancel = useCallback(() => {
    close()
  }, [close])

  return {
    // State
    visible,
    editing,
    submitting,
    form,

    // Actions
    open,
    close,
    submit,
    handleOk,
    handleCancel,
  }
}

export default useModal
